module KrangExt

using Synchray
using Synchray: GRLens, Ray, RayGR2, BLCoords, GeodesicRayState, CameraKerrGR,
	FourPosition, FourFrequency, FourVelocity,
	photon_k, _perpendicular_basis_vector, _captured_ray, camera_ray, direction3, frequency,
	_init_acc, _integrate_ray_step, _postprocess_acc, _rt_step_scalar,
	AccValue, Intensity, IntensityIQU, OpticalDepth, SpectralIndex,
	AbstractMedium, CombinedMedium,
	@swiz
using Krang: Krang, Kerr, SlowLightIntensityPixel, emission_coordinates, p_bl_d, metric_uu
using LinearAlgebra: cross, dot, norm
using StaticArrays: SVector, SMatrix

"""
	compute_deflection_map(lens::GRLens, cam)

Compute a gravitational deflection map for each pixel of `cam`.

Returns a collection with the same topology as `cam.xys`, where each element is either:
- a `Ray` — the outgoing asymptotic ray after BH deflection
- `nothing` — the ray is captured by the BH

# Arguments
- `lens`: black-hole lens specification
- `cam`: a flat-space camera supporting `camera_ray(cam, uv)`
"""
function Synchray.compute_deflection_map(
	lens::Synchray.GRLens,
	cam::Union{Synchray.CameraOrtho, Synchray.CameraPerspective},
)
	metric = Kerr(lens.spin)

	cam.mapfunc(cam.xys) do uv
		_deflected_ray(metric, camera_ray(cam, uv), lens)
	end
end

function Synchray.RayGR2(ray_in::Synchray.Ray, lens::Synchray.GRLens)
	ray_out = _deflected_ray(Kerr(lens.spin), ray_in, lens)
	Synchray.RayGR2(ray_in, ray_out, lens.bh_position)
end

function _deflected_ray(metric, ray_in, lens)
	n_to_bh = direction3(ray_in)
	impact = _impact_vector(ray_in, lens.bh_position)
	αhat, βhat, θ_obs = _observer_screen_basis(n_to_bh)
	α = dot(impact, αhat) / lens.bh_rg
	β = dot(impact, βhat) / lens.bh_rg
	pix = _krang_pixel(metric, α, β, θ_obs)
	R = _krang_to_lab_rotation(αhat, βhat, n_to_bh, θ_obs)
	_outgoing_ray_from_pixel(pix, ray_in, R, lens.bh_position, lens.bh_rg, lens.τ_frac)
end

"""Perpendicular vector from the BH to the ray line (impact parameter vector)."""
@inline _impact_vector(ray_in, bh_position) = begin
	anchor = (@swiz ray_in.x0.xyz) - bh_position
	n = direction3(ray_in)
	anchor - dot(anchor, n) * n
end

"""Observer screen basis and inclination for Krang, from the photon direction toward the observer."""
function _observer_screen_basis(n)
	spin_axis = SVector(0.0, 0.0, 1.0)
	# αhat: horizontal screen axis, ⊥ both spin axis and viewing direction.
	# Fallback for face-on viewing (n ∥ spin axis) where cross product vanishes.
	α = cross(spin_axis, n)
	if norm(α) < 1e-8
		α_fallback = SVector(1.0, 0.0, 0.0)
		α = α_fallback - dot(α_fallback, n) * n
	end
	αhat = α / norm(α)
	# βhat: vertical screen axis, completes the (αhat, βhat, n) right-handed triad.
	βhat = cross(n, αhat)
	# θ_obs: observer inclination = angle between viewing direction and spin axis.
	# Clamped away from poles to avoid Krang singularities.
	θ_eps = 1e-6
	θ_obs = clamp(acos(clamp(dot(n, spin_axis), -1.0, 1.0)), θ_eps, π - θ_eps)
	αhat, βhat, θ_obs
end

function _krang_pixel(metric, α, β, θ_obs)
	α = Float64(α)
	β = Float64(β)
	# Krang has a sign(β)=0 degeneracy there that pins θ=π/2.
	β = iszero(β) ? 1e-8 : β
	SlowLightIntensityPixel(metric, α, β, θ_obs)
end

"""Rotation matrix from Krang's Cartesian frame to Synchray's lab frame.
Chains: Krang Cartesian →(M_krang⁻¹) screen coords →(M_lab) lab frame."""
function _krang_to_lab_rotation(αhat, βhat, n_to_bh, θ_obs)
	sθ, cθ = sincos(θ_obs)
	# Krang's screen basis in Krang Cartesian (observer in y-z plane at inclination θ):
	M_krang = hcat(
		SVector(1.0, 0.0, 0.0),   # α̂: Krang α maps to +x̂ in BL Cartesian
		SVector(0.0, cθ, sθ),     # β̂: vertical screen axis, tilted by inclination
		SVector(0.0, -sθ, cθ),    # n̂: BH-to-observer direction
	)
	# Same screen basis in lab frame (from _observer_screen_basis):
	M_lab = hcat(αhat, βhat, n_to_bh)
	# M_krang' = M_krang⁻¹ (orthogonal), so R maps Krang Cartesian → lab
	M_lab * M_krang'
end

"""Construct the outgoing flat-space ray from a Krang geodesic solution.
Samples two points on the outgoing asymptote to get direction, then transforms to lab frame."""
function _outgoing_ray_from_pixel(pix, ray_in, R, bh_position, bh_rg, τ_frac)
	_is_captured(pix) && return _captured_ray(ray_in)

	τ_total = pix.total_mino_time
	if isnan(τ_total) || τ_total ≤ 0
		return _captured_ray(ray_in)
	end

	# Sample two points on the outgoing geodesic to determine direction.
	# τ_frac selects points on the asymptote where the ray is already straight.
	τ1 = τ_frac[1] * τ_total
	τ2 = τ_frac[2] * τ_total

	t1, r1, θ1, φ1, _, _, ok1 = emission_coordinates(pix, τ1)
	t2, r2, θ2, φ2, _, _, ok2 = emission_coordinates(pix, τ2)

	if !ok1 || !ok2
		return _captured_ray(ray_in)
	end

	# BL → Cartesian in Krang frame
	x1_krang = _bl_to_cartesian(r1, θ1, φ1)
	x2_krang = _bl_to_cartesian(r2, θ2, φ2)

	# Photon travel direction: from source toward BH (x1 is closer to BH than x2)
	dir = x1_krang - x2_krang
	n = norm(dir)
	if n < 1e-10
		return _captured_ray(ray_in)
	end
	n_out_krang = dir / n

	# Transform to lab: rotate, scale by bh_rg (Krang uses dimensionless units), translate to BH position
	n_out = R * n_out_krang
	x1_lab = bh_position + bh_rg * (R * x1_krang)

	k_out = photon_k(frequency(ray_in), n_out)
	# XXX: arbitrary e1 and e2 for now, so polarization is wrong! need to properly transport basis on deflection
	e1_out = _perpendicular_basis_vector(n_out, ray_in.e2)
	e2_out = cross(n_out, e1_out)
	Ray(
		FourPosition(ray_in.x0.t, x1_lab), k_out,
		e1_out, e2_out,
		ray_in.nz, ray_in.light, ray_in.s_max)
end

"""Check if a photon is captured by the BH: complex roots in the radial potential mean no turning point exists."""
function _is_captured(pix)
	roots = pix.roots
	any(abs(imag(r)) > 1e-6 for r in roots)
end

"""Boyer-Lindquist (r, θ, φ) → Cartesian (x, y, z) in Krang's frame.
BL: r = radial, θ = polar from spin axis, φ = azimuthal.
Cartesian: z = spin axis, x-y = equatorial plane, x = φ=0."""
_bl_to_cartesian(r, θ, φ) = SVector(
	r * sin(θ) * cos(φ),
	r * sin(θ) * sin(φ),
	r * cos(θ),
)


# ============================================================================
# Proper GR geodesic integration
# ============================================================================

"""
    KerrGeodesicPath

Wraps a Krang SlowLightIntensityPixel with coordinate transformation data.
Provides methods to step along the geodesic and convert BL → lab Cartesian.
"""
struct KerrGeodesicPath{T, TPix<:SlowLightIntensityPixel}
	metric::Kerr{T}
	pixel::TPix
	η::T                                    # reduced Carter constant
	λ::T                                    # reduced angular momentum
	ν_obs::T                                # observer frequency
	total_mino_time::T
	R_krang_to_lab::SMatrix{3,3,Float64,9}  # rotation: Krang Cartesian → lab frame
	bh_position::SVector{3,Float64}
	bh_rg::Float64
end

function KerrGeodesicPath(metric::Kerr, pix::SlowLightIntensityPixel, ν_obs,
		R_krang_to_lab, bh_position, bh_rg)
	KerrGeodesicPath(metric, pix, pix.η, pix.λ, Float64(ν_obs),
		pix.total_mino_time, R_krang_to_lab, bh_position, bh_rg)
end

"""
Convert BLCoords to lab-frame FourPosition.

Krang works in geometrized units where M = 1, so all coordinates are dimensionless
multiples of M. The physical length scale is set by `bh_rg` (gravitational radius, GM/c²).

**Spatial**: BL (r,θ,φ) → quasi-Cartesian Kerr-Schild (accounts for the BL→KS azimuthal
correction `φ_KS = φ_BL - a/(2√(1-a²)) ln|(r-r₊)/(r-r₋)|`) → rotate to lab frame
→ scale by bh_rg → translate to BH position.  Using Kerr-Schild instead of naive
spherical conversion is essential for spin > 0 where BL φ is frame-dragged.

**Time**: Krang's `emission_time_regularized` is the BL coordinate time integral along the
geodesic with logarithmic divergences subtracted (see Gralla & Lupsasca 2020).  It is **not**
the absolute BL coordinate time, but the finite part that encodes time delays between
different emission points.  Multiplied by bh_rg to convert to lab time units.
For `FastLight` models the time component is ignored; for `SlowLight` time-dependent models
this regularized time correctly captures relative delays between emission events.
"""
@inline function four_position_lab(path::KerrGeodesicPath, bl::BLCoords)
	x_ks = Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(path.metric, bl.r, bl.θ, bl.φ)
	x_krang = SVector(x_ks...)
	FourPosition(bl.t * path.bh_rg, path.bh_position + path.bh_rg * (path.R_krang_to_lab * x_krang))
end

"""
Convert BL contravariant spatial momentum to a Cartesian direction in Krang frame.

Uses Kerr-Schild φ: the Cartesian coordinates are functions of (r, θ, φ_KS) where
φ_KS = φ_BL + correction(r). The chain rule gives an extra dr contribution to dx/dφ_BL
from dφ_KS/dr. This is accounted for by computing the effective p^φ_KS contribution
via `p_φ_eff = p_φ + dφ_KS/dr * p_r`, then using the standard spherical→Cartesian Jacobian
with φ_KS.
"""
@inline function _bl_upper_to_cartesian_direction(metric::Kerr, p_r, p_θ, p_φ, r, θ, φ)
	a = metric.spin
	# φ_KS = φ_BL - a/(2√(1-a²)) ln|(r-r₊)/(r-r₋)| + atan(a,r)
	# dφ_KS/dr = log derivative + atan derivative
	temp = √(1 - a^2)
	rp = 1 + temp
	rm = 1 - temp
	dlog_dr = 1/(r - rp) - 1/(r - rm)           # d/dr ln|(r-r₊)/(r-r₋)|
	datan_dr = -a / (r^2 + a^2)                  # d/dr atan(a, r)
	dφ_KS_dr = -a/(2*temp) * dlog_dr + datan_dr
	# Effective p^φ in KS frame: chain rule dx/dr includes dφ_KS/dr contribution
	p_φ_eff = p_φ + dφ_KS_dr * p_r

	# Use KS azimuthal angle for the Jacobian
	φ_KS = Krang.ϕ_kerr_schild(metric, r, φ)
	sθ, cθ = sincos(θ)
	sφ, cφ = sincos(φ_KS)
	# Standard spherical→Cartesian Jacobian with φ_KS:
	# dx/dr = sinθ cosφ_KS,  dx/dθ = r cosθ cosφ_KS,  dx/dφ_KS = -r sinθ sinφ_KS
	vx = p_r * sθ * cφ + p_θ * r * cθ * cφ - p_φ_eff * r * sθ * sφ
	vy = p_r * sθ * sφ + p_θ * r * cθ * sφ + p_φ_eff * r * sθ * cφ
	vz = p_r * cθ      - p_θ * r * sθ
	v = SVector(vx, vy, vz)
	v / norm(v)
end

"""
Build a GeodesicRayState at a point on the geodesic.
Computes the local photon FourFrequency in lab Cartesian from BL momentum.
"""
@inline function geodesic_ray_state(path::KerrGeodesicPath, bl::BLCoords)
	# 1. BL photon momentum covector
	p_d = p_bl_d(path.metric, bl.r, bl.θ, path.η, path.λ, bl.νr, bl.νθ)

	# 2. Raise indices: p^μ = g^{μν} p_ν
	g_uu = metric_uu(path.metric, bl.r, bl.θ)
	p_u = g_uu * p_d

	# 3. BL contravariant spatial → KS Cartesian direction → lab frame
	n_krang = _bl_upper_to_cartesian_direction(path.metric, p_u[2], p_u[3], p_u[4], bl.r, bl.θ, bl.φ)
	n_lab = path.R_krang_to_lab * n_krang

	# 4. FourFrequency in lab Cartesian, with correct GR local energy.
	#    p_u[1] = p^t is the local photon energy per unit energy-at-infinity.
	#    Equals 1 at r→∞, grows near BH (gravitational blueshift for infalling photon).
	#    Embedding it into k ensures that lorentz_unboost gives the correct comoving
	#    frequency including gravitational redshift, not just SR Doppler.
	k = photon_k(path.ν_obs * p_u[1], n_lab)

	# 5. Polarization basis: placeholder (proper Walker-Penrose transport is Phase 3)
	e1 = _perpendicular_basis_vector(n_lab, SVector(0.0, 1.0, 0.0))
	e2 = cross(n_lab, e1)

	GeodesicRayState(k, e1, e2)
end

"""
    integrate_geodesic(obj, path, nτ, τ_range, what)

Integrate radiative transfer along a Kerr geodesic.

Steps in Mino time τ, converts BL→Cartesian at each step, and calls
the standard `_integrate_ray_step` with a duck-typed `GeodesicRayState`.
"""
function integrate_geodesic(obj::AbstractMedium, path::KerrGeodesicPath, nτ::Int,
		τ_range::NTuple{2,Float64}, what=Synchray.Intensity())
	ν = path.ν_obs
	a = path.metric.spin
	τ_total = path.total_mino_time

	# Handle fully invalid geodesics (no usable data at all)
	if isnan(τ_total) || τ_total ≤ 0
		return _postprocess_acc(_init_acc(typeof(what), ν), ν, what)
	end

	# Captured rays still have valid geodesic data along the incoming path
	# (camera → horizon). Krang's total_mino_time covers the full infall;
	# emission_coordinates returns ok=false once the ray hits the horizon,
	# so we naturally stop accumulating there.
	τ_min = τ_range[1] * τ_total
	τ_max = τ_range[2] * τ_total
	τs = range(τ_min, τ_max, length=nτ)
	Δτ = step(τs)

	acc = _init_acc(typeof(what), ν)

	r_horizon = Krang.horizon(path.metric)

	for τ in τs
		t_em, r, θ, φ, νr, νθ, ok = emission_coordinates(path.pixel, τ)
		!ok && continue
		# Skip points at or below the horizon — physics is unobservable there
		r ≤ r_horizon && continue
		# Skip numerically degenerate points (NaN θ from Krang β≈0 degeneracy)
		isnan(θ) && continue
		bl = BLCoords(t_em, r, θ, φ, νr, νθ)

		# Convert to Cartesian FourPosition
		x4 = four_position_lab(path, bl)

		# Affine parameter step: dλ = Σ dτ, scaled by bh_rg for physical units
		Σ = bl.r^2 + a^2 * cos(bl.θ)^2
		Δλ = Σ * Δτ * path.bh_rg

		# Build local photon state
		state = geodesic_ray_state(path, bl)

		# Reuse existing RT step — medium sees FourPosition + duck-typed ray
		acc = _integrate_ray_step(acc, obj, x4, state, Δλ)
	end

	_postprocess_acc(acc, ν, what)
end

# CombinedMedium: evaluate all sub-media at each geodesic point
function integrate_geodesic(cm::CombinedMedium, path::KerrGeodesicPath, nτ::Int,
		τ_range::NTuple{2,Float64}, what=Synchray.Intensity())
	ν = path.ν_obs
	a = path.metric.spin
	τ_total = path.total_mino_time

	if isnan(τ_total) || τ_total ≤ 0
		return _postprocess_acc(_init_acc(typeof(what), ν), ν, what)
	end

	τ_min = τ_range[1] * τ_total
	τ_max = τ_range[2] * τ_total
	τs = range(τ_min, τ_max, length=nτ)
	Δτ = step(τs)

	acc = _init_acc(typeof(what), ν)
	r_horizon = Krang.horizon(path.metric)

	for τ in τs
		t_em, r, θ, φ, νr, νθ, ok = emission_coordinates(path.pixel, τ)
		!ok && continue
		r ≤ r_horizon && continue
		isnan(θ) && continue
		bl = BLCoords(t_em, r, θ, φ, νr, νθ)

		x4 = four_position_lab(path, bl)
		Σ = bl.r^2 + a^2 * cos(bl.θ)^2
		Δλ = Σ * Δτ * path.bh_rg
		state = geodesic_ray_state(path, bl)

		for obj in cm.objects
			acc = _integrate_ray_step(acc, obj, x4, state, Δλ)
		end
	end

	_postprocess_acc(acc, ν, what)
end


"""
    render(cam::CameraKerrGR, obj, what)

Render a scene by tracing Kerr geodesics through the medium.
Per-pixel: computes Krang screen coordinates and inclination from ray direction,
creates a SlowLightIntensityPixel, and integrates RT along the geodesic.
"""
Synchray.render(cam::CameraKerrGR, obj::AbstractMedium, what=Synchray.Intensity()) = let
	metric = Kerr(cam.metric_spin)
	(; camera, bh_position, bh_rg, nτ, τ_range) = cam

	camera.mapfunc(camera.xys) do uv
		ray_in = camera_ray(camera, uv)

		# Per-pixel: direction, impact parameter, screen basis, inclination
		n = direction3(ray_in)
		impact = _impact_vector(ray_in, bh_position)
		αhat, βhat, θ_obs = _observer_screen_basis(n)
		α = dot(impact, αhat) / bh_rg
		β = dot(impact, βhat) / bh_rg
		R = _krang_to_lab_rotation(αhat, βhat, n, θ_obs)

		# Create Krang pixel (caches elliptic integrals)
		pix = _krang_pixel(metric, α, β, θ_obs)

		# Build geodesic path
		path = KerrGeodesicPath(metric, pix, frequency(ray_in), R, bh_position, bh_rg)

		# Integrate RT along geodesic
		integrate_geodesic(obj, path, nτ, τ_range, what)
	end
end

end

module KrangExt

using Synchray
using Synchray: GRLens, Ray, RayGR2, FourPosition, photon_k, _perpendicular_basis_vector, _captured_ray, camera_ray, direction3, frequency, @swiz
using Krang: Krang, Kerr, SlowLightIntensityPixel, emission_coordinates
using LinearAlgebra: cross, dot, norm

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
		ray_in.nz, ray_in.light)
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

end

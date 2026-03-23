module KrangExt

using Synchray
using Synchray: Ray, FourPosition, photon_k, _perpendicular_basis_vector
using Krang: Krang, Kerr, SlowLightIntensityCamera, emission_coordinates

"""
	compute_deflection_map(spin, θ_obs, cam::Synchray.CameraOrtho; τ_frac=(0.97, 0.99))

Compute a gravitational deflection map for each pixel of `cam`.

Returns an array (same shape as `cam.xys`) where each element is either:
- a `Ray` — the outgoing asymptotic ray after BH deflection
- `nothing` — the ray is captured by the BH

# Arguments
- `spin`: BH spin parameter a/M (use ≈0 for Schwarzschild, e.g. `1e-10`)
- `θ_obs`: observer inclination angle (radians, 0 = face-on, π/2 = edge-on)
- `cam`: a `CameraOrtho` defining the pixel grid. Screen coordinates `(u, v)` map to
  Krang screen coordinates `(α, β)` via the camera's `(e1, e2)` basis.
- `τ_frac`: two fractions of `total_mino_time` at which to evaluate the outgoing
  asymptote (default `(0.97, 0.99)` — both at large r past the BH).
"""
function Synchray.compute_deflection_map(spin, θ_obs, cam::Synchray.CameraOrtho; τ_frac=(0.97, 0.99))
	metric = Kerr(Float64(spin))

	# Build Krang camera matching the pixel grid
	αs = map(uv -> Float64(uv[1]), cam.xys)
	βs = map(uv -> Float64(uv[2]), cam.xys)
	αmin, αmax = extrema(αs)
	βmin, βmax = extrema(βs)
	res = size(cam.xys, 1)

	krang_cam = SlowLightIntensityCamera(metric, Float64(θ_obs), αmin, αmax, βmin, βmax, res)

	map(krang_cam.screen.pixels) do pix
		_outgoing_ray_from_pixel(pix, cam, τ_frac)
	end
end

function _outgoing_ray_from_pixel(pix, cam, τ_frac)
	if _is_captured(pix)
		return nothing
	end

	τ_total = pix.total_mino_time
	if isnan(τ_total) || τ_total ≤ 0
		return nothing
	end

	# Evaluate geodesic at two points on the outgoing asymptote
	τ1 = τ_frac[1] * τ_total
	τ2 = τ_frac[2] * τ_total

	t1, r1, θ1, φ1, _, _, ok1 = emission_coordinates(pix, τ1)
	t2, r2, θ2, φ2, _, _, ok2 = emission_coordinates(pix, τ2)

	if !ok1 || !ok2
		return nothing
	end

	# BL → Cartesian
	x1 = _bl_to_cartesian(r1, θ1, φ1)
	x2 = _bl_to_cartesian(r2, θ2, φ2)

	dir = x2 - x1
	n = norm(dir)
	if n < 1e-10
		return nothing
	end
	n_out = SVector{3}(dir / n)

	# Construct a full Ray
	k_out = photon_k(cam.ν, n_out)
	e1_out = _perpendicular_basis_vector(n_out, cam.e2)
	e2_out = cross(n_out, e1_out)
	Ray(FourPosition(cam.t, x1), k_out, e1_out, e2_out, cam.nz, cam.light)
end

function _is_captured(pix)
	roots = pix.roots
	any(abs(imag(r)) > 1e-6 for r in roots)
end

_bl_to_cartesian(r, θ, φ) = SVector(
	r * sin(θ) * cos(φ),
	r * sin(θ) * sin(φ),
	r * cos(θ),
)

end

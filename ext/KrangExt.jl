module KrangExt

using Synchray
using Synchray: Ray, FourPosition, photon_k, _perpendicular_basis_vector
using Krang: Krang, Kerr, SlowLightIntensityCamera, emission_coordinates

"""
	compute_deflection_map(spin, θ_obs, cam::CameraOrtho; bh_position=zeros(3), bh_rg=1.0, τ_frac=(0.97, 0.99))

Compute a gravitational deflection map for each pixel of `cam`.

Returns an array (same shape as `cam.xys`) where each element is either:
- a `Ray` — the outgoing asymptotic ray after BH deflection
- `nothing` — the ray is captured by the BH

# Arguments
- `spin`: BH spin parameter a/M
- `θ_obs`: observer inclination angle (radians, 0 = face-on, π/2 = edge-on)
- `cam`: a `CameraOrtho` defining the pixel grid
- `bh_position`: BH spatial position in lab coordinates (default: origin)
- `bh_rg`: gravitational radius GM/c² (default: 1). Krang coordinates are in units of r_g,
  so outgoing ray positions are `bh_position + bh_rg * krang_xyz`.
- `τ_frac`: two fractions of `total_mino_time` for outgoing asymptote evaluation
"""
function Synchray.compute_deflection_map(spin, θ_obs, cam::Synchray.CameraOrtho;
		bh_position=zero(SVector{3,Float64}), bh_rg=1.0, τ_frac=(0.97, 0.99))
	metric = Kerr(Float64(spin))

	# Krang operates in units of r_g. Convert pixel offsets from lab to Krang units.
	αs = map(uv -> Float64(uv[1]) / bh_rg, cam.xys)
	βs = map(uv -> Float64(uv[2]) / bh_rg, cam.xys)
	αmin, αmax = extrema(αs)
	βmin, βmax = extrema(βs)
	res = size(cam.xys, 1)

	krang_cam = SlowLightIntensityCamera(metric, Float64(θ_obs), αmin, αmax, βmin, βmax, res)

	map(krang_cam.screen.pixels) do pix
		_outgoing_ray_from_pixel(pix, cam, bh_position, bh_rg, τ_frac)
	end
end

function _outgoing_ray_from_pixel(pix, cam, bh_position, bh_rg, τ_frac)
	if _is_captured(pix)
		return nothing
	end

	τ_total = pix.total_mino_time
	if isnan(τ_total) || τ_total ≤ 0
		return nothing
	end

	τ1 = τ_frac[1] * τ_total
	τ2 = τ_frac[2] * τ_total

	t1, r1, θ1, φ1, _, _, ok1 = emission_coordinates(pix, τ1)
	t2, r2, θ2, φ2, _, _, ok2 = emission_coordinates(pix, τ2)

	if !ok1 || !ok2
		return nothing
	end

	# BL → Cartesian in Krang units, then scale to lab coordinates
	x1_krang = _bl_to_cartesian(r1, θ1, φ1)
	x2_krang = _bl_to_cartesian(r2, θ2, φ2)

	# Direction is scale-invariant
	dir = x2_krang - x1_krang
	n = norm(dir)
	if n < 1e-10
		return nothing
	end
	n_out = SVector{3}(dir / n)

	# Anchor in lab coordinates
	x1_lab = bh_position + bh_rg * x1_krang

	k_out = photon_k(cam.ν, n_out)
	e1_out = _perpendicular_basis_vector(n_out, cam.e2)
	e2_out = cross(n_out, e1_out)
	Ray(FourPosition(cam.t, x1_lab), k_out, e1_out, e2_out, cam.nz, cam.light)
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

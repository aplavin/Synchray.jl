module KrangExt

using Synchray
using Synchray: Ray, FourPosition, photon_k, _perpendicular_basis_vector
using Krang: Krang, Kerr, SlowLightIntensityPixel, emission_coordinates
using LinearAlgebra: normalize, cross, norm

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

function _krang_pixel(metric, uv, bh_rg, θ_obs)
	α = Float64(uv[1]) / bh_rg
	β = Float64(uv[2]) / bh_rg
	# Krang has a sign(β)=0 degeneracy there that pins θ=π/2.
	β = iszero(β) ? 1e-8 : β
	SlowLightIntensityPixel(metric, α, β, θ_obs)
end

"""Rotation matrix from Krang's Cartesian frame to Synchray's lab frame."""
function _krang_to_lab_rotation(θ_obs, cam)
	# Krang screen basis in Krang Cartesian:
	#   LOS (toward BH) = (0, sinθ, -cosθ)
	#   α̂ (horizontal)  = (-1, 0, 0)
	#   β̂ (vertical)    = (0, cosθ, sinθ)
	sθ, cθ = sincos(θ_obs)
	LOS = SVector(0.0, sθ, -cθ)
	α̂ = SVector(-1.0, 0.0, 0.0)
	β̂ = SVector(0.0, cθ, sθ)

	# Synchray camera basis: (cam.e1, cam.e2, cam.n)
	# R maps Krang basis → Synchray basis:
	#   R * α̂ = cam.e1,  R * β̂ = cam.e2,  R * (-LOS) = cam.n
	M_krang = hcat(α̂, β̂, -LOS)
	M_cam = hcat(cam.e1, cam.e2, cam.n)
	# Both are orthonormal, so R = M_cam * M_krang' (= M_cam * M_krang⁻¹)
	M_cam * M_krang'
end

function _outgoing_ray_from_pixel(pix, ray_in, R, bh_position, bh_rg, τ_frac)
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

	# BL → Cartesian in Krang frame
	x1_krang = _bl_to_cartesian(r1, θ1, φ1)
	x2_krang = _bl_to_cartesian(r2, θ2, φ2)

	# Direction in Krang frame (scale-invariant)
	dir = x2_krang - x1_krang
	n = norm(dir)
	if n < 1e-10
		return nothing
	end
	n_out_krang = SVector{3}(dir / n)

	# Rotate from Krang frame to Synchray lab frame
	n_out = SVector{3}(R * n_out_krang)
	x1_lab = bh_position + bh_rg * SVector{3}(R * x1_krang)

	k_out = photon_k(frequency(ray_in), n_out)
	e1_out = _perpendicular_basis_vector(n_out, ray_in.e2)
	e2_out = cross(n_out, e1_out)
	Ray(FourPosition(ray_in.x0.t, x1_lab), k_out, e1_out, e2_out, ray_in.nz, ray_in.light)
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

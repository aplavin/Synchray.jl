struct AngleTrigCached{T}
	tan::T
	cos::T
end
Base.tan(x::AngleTrigCached) = x.tan
Base.cos(x::AngleTrigCached) = x.cos
AngleTrigCached_fromangle(φ) = AngleTrigCached(tan(φ), cos(φ))

struct PowerLawS{Texp,Tval,Ts0}
	exp::Texp
	val0::Tval
	s0::Ts0
end

PowerLawS(exp; val0, s0=one(val0)) = PowerLawS(exp, val0, s0)

@inline (pl::PowerLawS)(s) = s > 0 ? pl.val0 * (s / pl.s0)^pl.exp : zero(float(pl.val0))

abstract type AbstractJetFieldDirection end

"""Marker direction model meaning: no ordered direction; `BFieldSpec.wrap` receives a scalar strength."""
struct ScalarField <: AbstractJetFieldDirection end

struct PoloidalField <: AbstractJetFieldDirection end
struct ToroidalField <: AbstractJetFieldDirection end
struct HelicalField{T} <: AbstractJetFieldDirection
	ψ::T
end

struct BFieldSpec{Tscale,Tdir,Twrap}
	scale::Tscale
	dir::Tdir
	wrap::Twrap
end

@inline jet_field_direction(::ScalarField, jet, x4) = 1

@inline jet_field_direction(::PoloidalField, jet, x4) = normalize(jet.axis)

"""
	jet_cylindrical_coords(axis, x4) -> (r, s, r_perp, ρ)

Compute cylindrical coordinates relative to jet axis.
Returns: spatial position `r`, axial distance `s`, perpendicular component `r_perp`, cylindrical radius `ρ`.
"""
@inline jet_cylindrical_coords(axis, x4) = begin
	r = @swiz x4.xyz
	s = dot(axis, r)
	r_perp = r - s * axis
	ρ = norm(r_perp)
	return (; r, s, r_perp, ρ)
end

@inline jet_field_direction(::ToroidalField, jet, x4) = begin
	axis = jet.axis
	(; r_perp, ρ) = jet_cylindrical_coords(axis, x4)
	iszero(ρ) && return zero(axis)
	e_r = r_perp / ρ
	return cross(axis, e_r)
end

@inline jet_field_direction(h::HelicalField, jet, x4) = begin
	axis = jet.axis
	(; r_perp, ρ) = jet_cylindrical_coords(axis, x4)
	e_phi = iszero(ρ) ? zero(axis) : cross(axis, r_perp / ρ)
	sψ, cψ = sincos(h.ψ)
	v = cψ * axis + sψ * e_phi
	return iszero(v) ? zero(v) : v / norm(v)
end

_rayz_cone_z_interval(axis, φj, ray::RayZ, s) = let
	# RayZ is a line of sight parameterized by z: r(z) = r0 + z*e_z.
	# Here we assume the standard camera convention: ray.x0.z == 0 and the cone apex is at the origin.
	@assert iszero(ray.x0.z)
	@assert 0 ≤ leftendpoint(s) ≤ rightendpoint(s)
	c2 = cos(φj)^2
	r0 = @swiz ray.x0.xyz
	az = axis.z

	# Cone inequality: (a⋅r)^2 ≥ c^2 * (r⋅r).
	# With r(z)=r0+z*e_z and r0⋅e_z=0:
	# (α0 + az*z)^2 - c^2*(|r0|^2 + z^2) >= 0
	α0 = dot(axis, r0)
	r02 = dot(r0, r0)
	A = az^2 - c2
	B = 2 * α0 * az
	C = α0^2 - c2 * r02

	FT = eltype(ray.x0) |> float
	emptyseg = let
		zref = FT(ray.x0.z)
		zref .. (zref - eps(zref))
	end
	infseg = FT(-Inf) .. FT(Inf)

	# 1) Restrict to truncation interval in s: s(z) ∈ s.
	zs = if iszero(az)
		(α0 ∈ s) ? infseg : emptyseg
	else
		smin, smax = endpoints(s)
		z1 = (smin - α0) / az
		z2 = (smax - α0) / az
		(min(z1, z2) .. max(z1, z2))
	end

	isempty(zs) && return zs

	# 2) Intersect with the infinite cone (already symmetric); half-cone handled above.
	# Solve A z^2 + B z + C ≥ 0 and clip to `zs`.
	z_cone = let
		D = B^2 - 4 * A * C
		if iszero(A)
			# Linear case: B z + C ≥ 0.
			if iszero(B)
				(C >= 0) ? infseg : emptyseg
			elseif B > 0
				(-C / B) .. infseg.right
			else
				infseg.left .. (-C / B)
			end
		elseif D < 0
			# No real roots: sign is constant (A and C have same sign when D<0).
			(C >= 0) ? infseg : emptyseg
		else
			z1 = (-B - √D) / (2 * A)
			z2 = (-B + √D) / (2 * A)
			zlo, zhi = minmax(z1, z2)
			if A < 0
				zlo .. zhi
			else
				# cone interior outside the roots
				@assert zlo ≤ 0 ≤ zhi  (zlo, zhi)
				left = infseg.left .. zlo
				right = zhi .. infseg.right
				# Prefer the side that overlaps `zs`
				L = zs ∩ left
				R = zs ∩ right
				@assert isempty(L) || isempty(R)
				!isempty(L) ? L : R
			end
		end
	end

	return zs ∩ z_cone
end

@inline _rotate_minimal(a::SVector{3}, b::SVector{3}, v::SVector{3}) = begin
	ET = float(promote_type(eltype(a), eltype(b), eltype(v)))

	# Assumes a and b are already unit vectors.
	c = dot(a, b)
	vx = cross(a, b)
	s = norm(vx)

	# Rodrigues rotation formula for the minimal rotation mapping â -> b̂.
	# Special-case (anti)parallel vectors for numerical stability.
	if s < √eps(s)
		if c > 0
			return ET.(v)
		else
			# 180° rotation: axis is not unique; pick a deterministic one perpendicular to â.
			# Since in this codepath â is typically ẑ, choosing ŷ keeps the "tilt about y" intuition.
			axis = abs(a.y) < 0.9 ? SVector{3,ET}(0, 1, 0) : SVector{3,ET}(1, 0, 0)
			k = normalize(cross(a, axis))
			# For θ = π: R(v) = v - 2 k×(k×v) = 2(k⋅v)k - v.
			return 2 * dot(k, v) * k - v
		end
	end

	k = vx / s
	# v_rot = v cosθ + (k×v) sinθ + k(k⋅v)(1-cosθ)
	return v * c + cross(k, v) * s + k * dot(k, v) * (1 - c)
end

"""
	jet_rotation_matrix(axis) -> SMatrix{3,3}

Return the 3×3 rotation matrix `R` implementing the jet coordinate convention:

- minimal rotation from the lab frame that maps the lab `ẑ` axis to `axis`.

The matrix columns are the lab-frame images of the jet basis vectors:

- first column: `ex` in lab coordinates
- second column: `ey` in lab coordinates
- third column: `ez == axis` in lab coordinates

Therefore:

- lab → jet coordinates: `r_jet = R' * r_lab`
- jet → lab coordinates: `r_lab = R * r_jet`
"""
@inline jet_rotation_matrix(axis::SVector{3}) = begin
	ẑ = SVector(0, 0, 1)
	x̂ = SVector(1, 0, 0)
	ŷ = SVector(0, 1, 0)

	ex = _rotate_minimal(ẑ, axis, x̂)
	ey = _rotate_minimal(ẑ, axis, ŷ)
	ez = axis

	return hcat(ex, ey, ez)
end

"""
	jet_basis(axis) -> (ex, ey, ez)

Construct an orthonormal right-handed spatial basis where `ez` points along `axis`.

Convention: **minimal rotation from lab frame**.
Apply the smallest 3D rotation that maps lab `ẑ` to the jet axis, and rotate lab `x̂`
and `ŷ` by the same rotation. This makes the jet frame as close to the lab frame as
possible while enforcing `ez ∥ axis`.
"""
@inline jet_basis(axis::SVector{3}) = eachcol(jet_rotation_matrix(axis))

"""
	lab_to_jet_coords(axis, r) -> SVector{3}

Convert a lab-frame spatial 3-position `r` to jet coordinates `(x, y, z)` defined by
`jet_basis(axis)`.
"""
@inline lab_to_jet_coords(axis::SVector{3}, r::SVector{3}) = jet_rotation_matrix(axis)' * r

"""
	jet_to_lab_coords(axis, rjet) -> SVector{3}

Convert a jet-frame spatial 3-position `rjet` back to lab coordinates.
This is the inverse of `lab_to_jet_coords`.
"""
@inline jet_to_lab_coords(axis::SVector{3}, rjet::SVector{3}) = jet_rotation_matrix(axis) * rjet

"""
    ray_in_jet_coords(ray, jet; z_range) -> Vector{SVector{3}}

Project a ray into jet coordinates.

Returns two 3-vectors in jet frame defining the line segment for the given `z_range` interval.
Use `@swiz` at call sites to extract desired components, e.g., `@swiz r.zx` for (s, h) coordinates.
"""
function ray_in_jet_coords(ray::RayZ, jet; z_range)
	map(endpoints(z_range)) do z
		r_lab = SVector(ray.x0.x, ray.x0.y, z)
		lab_to_jet_coords(jet, r_lab)
	end |> collect |> StructArray
end

"""
    camera_band_in_jet_plane(cam, jet; z_range) -> Vector{SVector{2}}

Project the camera's field of view into jet-plane coordinates (s, h).

Returns four corners of the band (polygon) swept by camera rays in the y=0 plane.

!!! note
    This function only uses the x-extent of the camera and assumes y=0, which is appropriate
    for 2D visualizations in the plane containing the jet axis and lab x-axis. It may not
    correctly represent the full camera FOV when the jet axis has significant y-component.
"""
@unstable function camera_band_in_jet_plane(cam::CameraZ, jet; z_range)
    xmin, xmax = extrema(xy -> xy.x, cam.xys)
    z1, z2 = endpoints(z_range)

    corners_lab = [
        SVector(xmin, 0.0, z1),
        SVector(xmax, 0.0, z1),
        SVector(xmax, 0.0, z2),
        SVector(xmin, 0.0, z2),
    ]

    map(corners_lab) do r_lab
        r_jet = lab_to_jet_coords(jet, r_lab)
        SVector(r_jet[3], r_jet[1])
    end
end


@kwdef struct ConicalJet{Ta,Tφ,Ts,Tne,TB,Tu,Tmodel} <: AbstractSynchrotronMedium
	axis::Ta
	φj::Tφ
	s::Ts
	ne::Tne
	B::TB
	speed_profile::Tu
	model::Tmodel
end

@unstable prepare_for_computations(obj::ConicalJet) = @p let
	obj
	@modify(AngleTrigCached_fromangle, __.φj)
	@modify(prepare_for_computations, __.model)
end

z_interval(obj::ConicalJet, ray::RayZ) = _rayz_cone_z_interval(obj.axis, obj.φj, ray, obj.s)

@inline is_inside_jet(jet::ConicalJet, x4::FourPosition) = begin
	(; s, ρ) = jet_cylindrical_coords(jet.axis, x4)
	return (s ∈ jet.s) && ρ ≤ s * tan(jet.φj)
end

jet_rotation_matrix(jet::ConicalJet) = jet_rotation_matrix(jet.axis)
jet_basis(jet::ConicalJet) = jet_basis(jet.axis)
lab_to_jet_coords(jet::ConicalJet, r::SVector{3}) = lab_to_jet_coords(jet.axis, r)
jet_to_lab_coords(jet::ConicalJet, rjet::SVector{3}) = jet_to_lab_coords(jet.axis, rjet)

@inline four_velocity(obj::ConicalJet, x4) = begin
	(; r, s, ρ) = jet_cylindrical_coords(obj.axis, x4)

	vhat = iszero(r) ? zero(r) : r / √(dot(r, r))
	tanφ = tan(obj.φj)
	η = ρ / (s * tanφ)

	(whichspeed, speed) = obj.speed_profile(η)
	return construct(FourVelocity, whichspeed => speed, direction => vhat)
end

@inline electron_density(obj::ConicalJet, x4) = begin
	s = dot(obj.axis, @swiz x4.xyz)
	return obj.ne(s)
end

@inline magnetic_field(obj::ConicalJet, x4) = begin
	s = dot(obj.axis, @swiz x4.xyz)
	b = obj.B.scale(s)
	b̂ = jet_field_direction(obj.B.dir, obj, x4)
	return obj.B.wrap(b * b̂)
end

@inline synchrotron_model(obj::ConicalJet) = obj.model

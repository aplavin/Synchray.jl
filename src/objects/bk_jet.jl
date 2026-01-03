struct AngleTrigCached{T}
	tan::T
	cos::T
end
Base.tan(x::AngleTrigCached) = x.tan
Base.cos(x::AngleTrigCached) = x.cos
AngleTrigCached_fromangle(φ) = AngleTrigCached(tan(φ), cos(φ))

@kwdef struct ConicalBKJet{Ta,Tφ,Ts,Ts0,Tne0,TB0,Tneexp,TBexp,Tu,Tmodel} <: AbstractSynchrotronMedium
	axis::Ta
	φj::Tφ
	s::Ts
	s0::Ts0
	ne0::Tne0
	B0::TB0
	ne_exp::Tneexp = -2
	B_exp::TBexp = -1
	speed_profile::Tu
	model::Tmodel = AngleAveragedPowerLawElectrons(p=2.5)
end

@unstable prepare_for_computations(obj::ConicalBKJet) = @p let
	obj
	@modify(AngleTrigCached_fromangle, __.φj)
	@modify(prepare_for_computations, __.model)
	modify(FixedExponent, __, @o _.ne_exp _.B_exp)
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

z_interval(obj::ConicalBKJet, ray::RayZ) = _rayz_cone_z_interval(obj.axis, obj.φj, ray, obj.s)

@inline _rotate_minimal(a::SVector{3}, b::SVector{3}, v::SVector{3}) = begin
	# Assumes a and b are already unit vectors.
	c = dot(a, b)
	vx = cross(a, b)
	s = norm(vx)

	# Rodrigues rotation formula for the minimal rotation mapping â -> b̂.
	# Special-case (anti)parallel vectors for numerical stability.
	if s < √eps(s)
		if c > 0
			return float(v)
		else
			# 180° rotation: axis is not unique; pick a deterministic one perpendicular to â.
			# Since in this codepath â is typically ẑ, choosing ŷ keeps the "tilt about y" intuition.
			axis = abs(a.y) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
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

- minimal rotation from the lab frame that maps the lab `ẑ` axis to `axis`.

The matrix columns are the lab-frame images of the jet basis vectors:

- first column: `ex` in lab coordinates
- second column: `ey` in lab coordinates
- third column: `ez == axis` in lab coordinates

Therefore:

- lab → jet coordinates: `r_jet = R' * r_lab`
- jet → lab coordinates: `r_lab = R * r_jet`
"""
@inline jet_rotation_matrix(axis::SVector{3}) = begin
	ẑ = SVector(0, 0, 1)
	x̂ = SVector(1, 0, 0)
	ŷ = SVector(0, 1, 0)

	ex = _rotate_minimal(ẑ, axis, x̂)
	ey = _rotate_minimal(ẑ, axis, ŷ)
	ez = axis
	
	return hcat(ex, ey, ez)
end
jet_rotation_matrix(jet::ConicalBKJet) = jet_rotation_matrix(jet.axis)

"""
	jet_basis(axis) -> (ex, ey, ez)

Construct an orthonormal right-handed spatial basis where `ez` points along `axis`.

Convention: **minimal rotation from lab frame**.
Apply the smallest 3D rotation that maps lab `ẑ` to the jet axis, and rotate lab `x̂`
and `ŷ` by the same rotation. This makes the jet frame as close to the lab frame as
possible while enforcing `ez ∥ axis`.
"""
@inline jet_basis(axis::SVector{3}) = eachcol(jet_rotation_matrix(axis))
jet_basis(jet::ConicalBKJet) = jet_basis(jet.axis)

"""
	lab_to_jet_coords(axis, r) -> SVector{3}

Convert a lab-frame spatial 3-position `r` to jet coordinates `(x, y, z)` defined by
`jet_basis(axis)`.
"""
@inline lab_to_jet_coords(axis::SVector{3}, r::SVector{3}) = jet_rotation_matrix(axis)' * r
lab_to_jet_coords(jet::ConicalBKJet, r::SVector{3}) = lab_to_jet_coords(jet.axis, r)

"""
	jet_to_lab_coords(axis, rjet) -> SVector{3}

Convert a jet-frame spatial 3-position `rjet` back to lab coordinates.
This is the inverse of `lab_to_jet_coords`.
"""
@inline jet_to_lab_coords(axis::SVector{3}, rjet::SVector{3}) = jet_rotation_matrix(axis) * rjet
jet_to_lab_coords(jet::ConicalBKJet, rjet::SVector{3}) = jet_to_lab_coords(jet.axis, rjet)

"""
	is_inside_jet(jet::ConicalBKJet, x4::FourPosition) -> Bool

Return whether the spacetime event `x4` lies inside the truncated conical jet volume.

The check uses only the spatial part `r = x4.xyz`:

- longitudinal coordinate `s = jet.axis ⋅ r` must lie in `jet.s` (truncation)
- transverse distance from the axis must satisfy `ρ ≤ s * tan(φj)` (cone boundary)

Boundary points count as inside.
"""
@inline is_inside_jet(jet::ConicalBKJet, x4::FourPosition) = begin
	r = @swiz x4.xyz
	s = dot(jet.axis, r)
	(s ∈ jet.s) || return false
	ρ = norm(r - s * jet.axis)
	return ρ ≤ s * tan(jet.φj)
end

@inline four_velocity(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)

	vhat = iszero(r) ? zero(r) : r / √(dot(r, r))
	# Transverse coordinate η in [0,1]
	ρ = norm(r - s * obj.axis)
	tanφ = tan(obj.φj)
	η = ρ / (s * tanφ)

	(whichspeed, speed) = obj.speed_profile(η)
	return construct(FourVelocity, whichspeed => speed, direction => vhat)
end

@inline electron_density(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	s > 0 ? obj.ne0 * (s / obj.s0)^(obj.ne_exp) : zero(float(obj.ne0))
end

@inline magnetic_field_strength(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	s > 0 ? obj.B0 * (s / obj.s0)^(obj.B_exp) : zero(float(obj.B0))
end

@inline synchrotron_model(obj::ConicalBKJet) = obj.model






"""
    ConicalBKJetWithPatterns

Wrapper medium that represents:

    (ConicalBKJet baseline) × (multiplicative pattern factors)

Design intent:

- Delegate geometry and kinematics to `base` (same `z_interval`, same `four_velocity`).
- Apply patterns only to comoving `electron_density` and `magnetic_field_strength`.
"""
struct ConicalBKJetWithPatterns{Tbase, Tpatterns} <: AbstractSynchrotronMedium
	base::Tbase
	patterns::Tpatterns
end

function ConicalBKJetWithPatterns(base::ConicalBKJet, patterns)
	for p in patterns
		validate_pattern(p, base)
	end
	return ConicalBKJetWithPatterns{typeof(base), typeof(patterns)}(base, patterns)
end

@unstable prepare_for_computations(obj::ConicalBKJetWithPatterns) = @p let
	obj
	@modify(prepare_for_computations, __.base)
end

@inline z_interval(obj::ConicalBKJetWithPatterns, ray::RayZ) = z_interval(obj.base, ray)
@inline four_velocity(obj::ConicalBKJetWithPatterns, x4) = four_velocity(obj.base, x4)

@inline is_inside_jet(obj::ConicalBKJetWithPatterns, x4::FourPosition) = is_inside_jet(obj.base, x4)
@inline jet_basis(obj::ConicalBKJetWithPatterns) = jet_basis(obj.base)
@inline lab_to_jet_coords(obj::ConicalBKJetWithPatterns, r::SVector{3}) = lab_to_jet_coords(obj.base, r)
@inline jet_to_lab_coords(obj::ConicalBKJetWithPatterns, rjet::SVector{3}) = jet_to_lab_coords(obj.base, rjet)

@inline electron_density(obj::ConicalBKJetWithPatterns, x4) =
	electron_density(obj.base, x4) * pattern_factor_ne(obj.patterns, x4, obj.base)

@inline magnetic_field_strength(obj::ConicalBKJetWithPatterns, x4) =
	magnetic_field_strength(obj.base, x4) * pattern_factor_B(obj.patterns, x4, obj.base)

@inline synchrotron_model(obj::ConicalBKJetWithPatterns) = synchrotron_model(obj.base)

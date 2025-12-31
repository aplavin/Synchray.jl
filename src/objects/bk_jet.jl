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
	model::Tmodel = PowerLawElectrons(p=2.5)
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

@inline four_velocity(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)

	vhat = r / √(dot(r, r))
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

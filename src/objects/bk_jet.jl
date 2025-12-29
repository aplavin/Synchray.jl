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

_inside_cone(axis::SVector{3}, φj, r::SVector{3}) = dot(axis, r)^2 >= cos(φj)^2 * dot(r, r)

_z_interval_from_s(axis::SVector{3}, s_interval, r0::SVector{3}) = begin
	α = dot(axis, r0)
	β = axis.z
	@assert !iszero(axis.z)
	zmin_s = (leftendpoint(s_interval) - α) / β
	zmax_s = (rightendpoint(s_interval) - α) / β
	return min(zmin_s, zmax_s) .. max(zmin_s, zmax_s)
end

z_interval(obj::ConicalBKJet, ray::RayZ) = begin
	φj = obj.φj

	r0 = @swizzle ray.x0.xyz
	z0 = ray.x0.z

	# Enforce one-nappe cone from apex: s >= 0 and within obj.s
	z_dom = _z_interval_from_s(obj.axis, obj.s, r0)
	isempty(z_dom) && return z_dom

	# Cone inequality along the ray:
	#   f(z) = (a⋅r(z))^2 - cos(φj)^2 * |r(z)|^2 ≥ 0,
	# with r(z) = r0 + z*(0,0,1).
	# Let α = a⋅r0 and β = a⋅e_z (= a_z). Then
	#   a⋅r(z) = α + β z
	#   |r(z)|^2 = |r0|^2 + 2 z0 z + z^2
	# so f(z) is a quadratic A z^2 + B z + C.
	α = dot(obj.axis, r0)
	β = obj.axis.z
	r02 = dot(r0, r0)
	c2 = cos(φj)^2
	A = β * β - c2
	B = 2 * (α * β - c2 * z0)
	C = α * α - c2 * r02

	T = promote_type(eltype(endpoints(z_dom)), typeof(A), typeof(B), typeof(C))
	Az = T(A)
	Bz = T(B)
	Cz = T(C)

	roots = T[]
	if iszero(Az)
		if !iszero(Bz)
			push!(roots, -Cz / Bz)
		end
	else
		D = Bz^2 - 4 * Az * Cz
		if D >= 0
			push!(roots, (-Bz - √D) / (2 * Az))
			push!(roots, (-Bz + √D) / (2 * Az))
		end
	end

	cut_vals = collect(endpoints(z_dom))
	for r in roots
		r ∈ z_dom && push!(cut_vals, r)
	end
	cut_vals = sort(unique(cut_vals))

	f(z) = (Az * z + Bz) * z + Cz
	s_at(z) = α + β * z

	inside_segments = ClosedInterval{T}[]
	for i in 1:(length(cut_vals) - 1)
		zab = cut_vals[i]..cut_vals[i + 1]
		zmid = mean(zab)
		if (f(zmid) >= 0) && (s_at(zmid) >= 0)
			push!(inside_segments, zab)
		end
	end

    @info "" z_dom repr(cut_vals) repr(roots) repr(inside_segments)
	if isempty(inside_segments)
		zz = float(z0)
		return zz .. (zz - eps(zz))
	end
	# Intersection should be connected for a single-nappe cone; merge defensively.
	z1 = minimum(leftendpoint.(inside_segments))
	z2 = maximum(rightendpoint.(inside_segments))
	return z1 .. z2
end

four_velocity(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	@assert s ≤ 0

	r2 = dot(r, r)
	vhat = r / sqrt(r2)

	# Transverse coordinate η in [0,1]
	ρ = norm(r - s * obj.axis)
	tanφ = tan(obj.φj)
	η = ρ / (s * tanφ)

	β = obj.speed_profile(η)
	return FourVelocity(β * vhat)
end

electron_density(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	(s ∈ obj.s) || return zero(s)
	_inside_cone(obj.axis, obj.φj, r) || return zero(s)
	s > 0 || return zero(s)
	return obj.ne0 * (s / obj.s0)^(obj.ne_exp)
end

magnetic_field_strength(obj::ConicalBKJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	(s ∈ obj.s) || return zero(s)
	_inside_cone(obj.axis, obj.φj, r) || return zero(s)
	s > 0 || return zero(s)
	return obj.B0 * (s / obj.s0)^(obj.B_exp)
end

synchrotron_model(obj::ConicalBKJet) = obj.model

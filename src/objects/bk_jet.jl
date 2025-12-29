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
	speed_kind::Symbol = :beta  # :beta or :gamma
	model::Tmodel = PowerLawElectrons(p=2.5)
end

_inside_cone(axis::SVector{3}, φj, r::SVector{3}) = begin
	s = dot(axis, r)
	s <= 0 && return false
	c = cos(φj)
	c2 = c * c
	return (s * s) >= c2 * dot(r, r)
end

_cone_boundary_poly(axis::SVector{3}, φj, x, y) = begin
	ax, ay, az = axis
	c = cos(φj)
	c2 = c * c
	s0 = ax * x + ay * y
	rxy2 = x * x + y * y
	A = az * az - c2
	B = 2 * s0 * az
	C = s0 * s0 - c2 * rxy2
	return A, B, C
end

_z_interval_from_s(axis::SVector{3}, s_interval, x, y, z0) = begin
	ax, ay, az = axis
	s0 = ax * x + ay * y + az * z0
	if iszero(az)
		return (s0 ∈ s_interval) ? (-Inf .. Inf) : (zero(z0) .. (zero(z0) - eps(float(z0))))
	end
	zmin_s = (leftendpoint(s_interval) - s0) / az
	zmax_s = (rightendpoint(s_interval) - s0) / az
	return min(zmin_s, zmax_s) .. max(zmin_s, zmax_s)
end

z_interval(obj::ConicalBKJet, ray::RayZ) = begin
	φj = obj.φj

	x = ray.x0.x
	y = ray.x0.y
	z0 = ray.x0.z

	# Enforce one-nappe cone from apex: s >= 0 and within obj.s
	z_s = _z_interval_from_s(obj.axis, obj.s, x, y, z0)
	if isempty(z_s)
		return z_s
	end
	# additionally enforce s >= 0 even if obj.s includes negatives
	z_spos = _z_interval_from_s(obj.axis, 0 .. Inf, x, y, z0)
	z_dom = intersect(z_s, z_spos)
	if isempty(z_dom)
		return z_dom
	end

	# Solve cone inequality in terms of absolute spatial z = z0 + z
	A, B, C = _cone_boundary_poly(obj.axis, φj, x, y)
	# Shift polynomial from Z (absolute z) to z (ray parameter): Z = z0 + z
	# f(Z) = A Z^2 + B Z + C; f(z0 + z) expands to:
	A′ = A
	B′ = 2 * A * z0 + B
	C′ = A * z0 * z0 + B * z0 + C

	zmin = leftendpoint(z_dom)
	zmax = rightendpoint(z_dom)

	T = promote_type(typeof(zmin), typeof(zmax), typeof(A′), typeof(B′), typeof(C′))
	Az = T(A′)
	Bz = T(B′)
	Cz = T(C′)

	cut_vals = T[T(zmin), T(zmax)]

	roots = T[]
	if iszero(Az)
		if !iszero(Bz)
			push!(roots, -Cz / Bz)
		end
	else
		D = Bz * Bz - 4 * Az * Cz
		if D >= 0
			√D = sqrt(D)
			push!(roots, (-Bz - √D) / (2 * Az))
			push!(roots, (-Bz + √D) / (2 * Az))
		end
	end

	for r in roots
		(r > zmin && r < zmax) && push!(cut_vals, r)
	end
	cut_vals = sort(unique(cut_vals))

	f(z) = (Az * z + Bz) * z + Cz
	s_at(z) = dot(obj.axis, SVector(x, y, z0 + z))

	inside_segments = Tuple{T, T}[]
	for i in 1:(length(cut_vals) - 1)
		za = cut_vals[i]
		zb = cut_vals[i + 1]
		zmid = (za + zb) / 2
		if (f(zmid) >= 0) && (s_at(zmid) >= 0)
			push!(inside_segments, (za, zb))
		end
	end

	if isempty(inside_segments)
		zz = float(z0)
		return zz .. (zz - eps(zz))
	end
	# Intersection should be connected for a single-nappe cone; merge defensively.
	z1 = minimum(first.(inside_segments))
	z2 = maximum(last.(inside_segments))
	return z1 .. z2
end

four_velocity(obj::ConicalBKJet, x4) = begin
	r = SVector(x4.x, x4.y, x4.z)
	s = dot(obj.axis, r)
	if !(s > 0)
		return FourVelocity(SVector(zero(s), zero(s), zero(s)))
	end

	r2 = dot(r, r)
	if iszero(r2)
		return FourVelocity(SVector(zero(s), zero(s), zero(s)))
	end
	vhat = r / sqrt(r2)

	# Transverse coordinate η in [0,1]
	ρ = norm(r - s * obj.axis)
	tanφ = tan(obj.φj)
	η = (tanφ == 0) ? zero(ρ) : clamp(ρ / (s * tanφ), zero(ρ), one(ρ))

	v = obj.speed_profile(η)
	β = obj.speed_kind === :gamma ? sqrt(max(zero(v), one(v) - inv(v * v))) : v
	β = clamp(β, zero(β), one(β) - eps(float(one(β))))
	return FourVelocity(β * vhat)
end

electron_density(obj::ConicalBKJet, x4) = begin
	r = SVector(x4.x, x4.y, x4.z)
	s = dot(obj.axis, r)
	(s ∈ obj.s) || return zero(s)
	_inside_cone(obj.axis, obj.φj, r) || return zero(s)
	s > 0 || return zero(s)
	return obj.ne0 * (s / obj.s0)^(-obj.ne_exp)
end

magnetic_field_strength(obj::ConicalBKJet, x4) = begin
	r = SVector(x4.x, x4.y, x4.z)
	s = dot(obj.axis, r)
	(s ∈ obj.s) || return zero(s)
	_inside_cone(obj.axis, obj.φj, r) || return zero(s)
	s > 0 || return zero(s)
	return obj.B0 * (s / obj.s0)^(-obj.B_exp)
end

synchrotron_model(obj::ConicalBKJet) = obj.model

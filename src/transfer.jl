struct AccValue{WHAT, T}
	value::T
end
AccValue{WHAT}(value) where {WHAT} = AccValue{WHAT, typeof(value)}(value)

_init_acc(::Type{Intensity}) = AccValue{Intensity}(0)
_init_acc(::Type{OpticalDepth}) = AccValue{OpticalDepth}(0)
_init_acc(::Type{T}) where {T<:Tuple} = map(_init_acc, fieldtypes(T))

_postprocess_acc(Iinv::AccValue{Intensity}, ν) = Iinv.value * ν^3
_postprocess_acc(τ::AccValue{OpticalDepth}, ν) = τ.value
_postprocess_acc(acc::Tuple, ν) = map(a -> _postprocess_acc(a, ν), acc)


_integrate_ray_step(acc::AccValue{Intensity}, obj, x4, k, Δλ, what::Intensity) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u)
	Jinv = emissivity_invariant(obj, x4, ν)
	Ainv = absorption_invariant(obj, x4, ν)

	Iinv = acc.value
	Δτ = Ainv * Δλ
	Iinv = if abs(Δτ) < 1e-8
        Iinv + (Jinv - Ainv * Iinv) * Δλ
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return AccValue{Intensity}(Iinv)
end

_integrate_ray_step(acc::AccValue{OpticalDepth}, obj, x4, k, Δλ, what::OpticalDepth) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u)
	Ainv = absorption_invariant(obj, x4, ν)

	Δτ = Ainv * Δλ
	return AccValue{OpticalDepth}(acc.value + Δτ)
end

_integrate_ray_step(acc::Tuple, obj, x4, k, Δλ, what::Tuple{Intensity,OpticalDepth}) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u)
	Jinv = emissivity_invariant(obj, x4, ν)
	Ainv = absorption_invariant(obj, x4, ν)

	Iinv = acc[1].value
	Δτ = Ainv * Δλ
	Iinv = if abs(Δτ) < 1e-8
        Iinv + (Jinv - Ainv * Iinv) * Δλ
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return (AccValue{Intensity}(Iinv),
			AccValue{OpticalDepth}(acc[2].value + Δτ))
end


integrate_ray(obj::AbstractMedium, ray::RayZ, what=Intensity()) = begin
	seg = z_interval(obj, ray)

	k = ray.k
	kz = k.z
    @assert k == SVector(kz, 0, 0, kz)
    k1 = k / kz

	zs = range(seg, isempty(seg) ? 0 : ray.nz)
	Δz = step(zs)
	Δλ = Δz / kz

	acc = _integrate_ray_step(_init_acc(typeof(what)), obj, ray.x0 + first(zs) * k1, k, Δλ, what)
	for z in zs[2:end]
        x = ray.x0 + z * k1
		acc = _integrate_ray_step(acc, obj, x, k, Δλ, what)
	end

	ν = photon_frequency(ray.k)
	return _postprocess_acc(acc, ν)
end

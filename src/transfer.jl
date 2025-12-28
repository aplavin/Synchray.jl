struct AccValue{WHAT, T}
	value::T
end
AccValue{WHAT}(value) where {WHAT} = AccValue{WHAT, typeof(value)}(value)

_init_acc(::Type{T}, ν0) where {T<:Tuple} = map(a -> _init_acc(a, ν0), fieldtypes(T))
_init_acc(::Type{Intensity}, ν0) = AccValue{Intensity}(0)
_init_acc(::Type{OpticalDepth}, ν0) = AccValue{OpticalDepth}(0)
struct AccSpectralIndex{TI, TS, DT}
	Iinv::TI
	s::TS
	DT::Type{DT}  # for type stability; this is the reason for a custom struct instead of a NamedTuple
end
_init_acc(::Type{Tuple{Intensity,SpectralIndex}}, ν0) = _init_acc(SpectralIndex, ν0)
_init_acc(::Type{SpectralIndex}, ν0) = begin
	ν0f = float(ν0)
	DT = typeof(ForwardDiff.Tag(_integrate_ray_step, typeof(ν0f)))
	s = ForwardDiff.Dual{DT}(one(ν0f), inv(ν0f))
	Iinv = ForwardDiff.Dual{DT}(zero(ν0f), zero(ν0f))
	AccSpectralIndex(Iinv, s, DT)
end


_postprocess_acc(acc::Tuple, ν, what::Tuple) = map((a, w) -> _postprocess_acc(a, ν, w), acc, what)
_postprocess_acc(Iinv::AccValue{Intensity}, ν, what::Intensity) = Iinv.value * ν^3
_postprocess_acc(τ::AccValue{OpticalDepth}, ν, what::OpticalDepth) = τ.value
_postprocess_acc(acc::AccSpectralIndex, ν, what::SpectralIndex) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	return (ν / I0) * dI
end
_postprocess_acc(acc::AccSpectralIndex, ν, what::Tuple{Intensity,SpectralIndex}) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	return I0, (ν / I0) * dI
end


_integrate_ray_step(acc::AccValue{Intensity}, obj, x4, k, Δλ) = begin
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

_integrate_ray_step(acc::AccValue{OpticalDepth}, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u)
	Ainv = absorption_invariant(obj, x4, ν)

	Δτ = Ainv * Δλ
	return AccValue{OpticalDepth}(acc.value + Δτ)
end

_integrate_ray_step(acc::AccSpectralIndex, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u) * acc.s
	Δλ′ = Δλ / acc.s
	Jinv = emissivity_invariant(obj, x4, ν)
	Ainv = absorption_invariant(obj, x4, ν)

	Iinv = acc.Iinv
	Δτ = Ainv * Δλ′
	Δτ0 = ForwardDiff.value(Δτ)
	Iinv = if abs(Δτ0) < 1e-8
		Iinv + (Jinv - Ainv * Iinv) * Δλ′
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return AccSpectralIndex(Iinv, acc.s, acc.DT)
end

_integrate_ray_step(acc::Tuple{AccValue{Intensity}, AccValue{OpticalDepth}}, obj, x4, k, Δλ) = begin
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

	acc = _integrate_ray_step(_init_acc(typeof(what), photon_frequency(ray)), obj, ray.x0 + first(zs) * k1, k, Δλ)
	for z in zs[2:end]
        x = ray.x0 + z * k1
		acc = _integrate_ray_step(acc, obj, x, k, Δλ)
	end

	ν = photon_frequency(ray)
	return _postprocess_acc(acc, ν, what)
end

_step_invariant(Iinv, Jinv, Ainv, Δλ) = begin
	(Ainv == 0 || Δλ == 0) && return Iinv + Jinv * Δλ
	Δτ = Ainv * Δλ
	if abs(Δτ) < 1e-8
        Iinv + (Jinv - Ainv * Iinv) * Δλ
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
end

_integrate_ray_step(Iinv, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	ν = measured_frequency(k, u)
	Jinv = emissivity_invariant(obj, x4, ν)
	Ainv = absorption_invariant(obj, x4, ν)
	return _step_invariant(Iinv, Jinv, Ainv, Δλ)
end

integrate_ray(obj::AbstractMedium, ray::RayZ) = begin
	seg = z_interval(obj, ray)

	k = ray.k
	kz = k.z
    @assert k == SVector(kz, 0, 0, kz)
    k1 = k / kz

	zs = range(seg, isempty(seg) ? 0 : ray.nz)
	Δz = step(zs)
	Δλ = Δz / kz

	Iinv = _integrate_ray_step(0, obj, ray.x0 + first(zs) * k1, k, Δλ)
	for z in zs[2:end]
        x = ray.x0 + z * k1
		Iinv = _integrate_ray_step(Iinv, obj, x, k, Δλ)
	end

	ν = photon_frequency(ray.k)
	return ν^3 * Iinv
end

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
	_step_invariant(Iinv, Jinv, Ainv, Δλ)
end

integrate_ray(obj::AbstractMedium, xy::SVector{2}; νcam, t0=0.0, nz=256) = begin
	seg = z_interval(obj, xy, νcam, t0)

    x0 = FourPosition(t0, xy..., 0)

	k = photon_k(νcam)
	kz = k.z
    @assert k == SVector(kz, 0, 0, kz)
    k1 = k / kz

	zs = range(seg, isempty(seg) ? 0 : nz)
	Δz = step(zs)
	Δλ = Δz / kz

	Iinv = _integrate_ray_step(0, obj, x0 + first(zs) * k1, k, Δλ)
	for z in zs[2:end]
        x = x0 + z * k1
		Iinv = _integrate_ray_step(Iinv, obj, x, k, Δλ)
	end

	νcam^3 * Iinv
end

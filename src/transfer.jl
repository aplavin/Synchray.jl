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

integrate_ray(obj::AbstractMedium, xy::SVector{2}; νcam, t0=0.0, nz=256) = begin
	seg = z_interval(obj, xy, νcam, t0)

	k = photon_k(νcam)
	kz = k.z
	Iinv = Iinv0

    zs = range(seg, isempty(seg) ? 0 : nz)
    Δz = step(zs)
	for z in zs
		t = t0 + z
		x4 = FourPosition(t, xy..., z)
		u = four_velocity(obj, x4)
		ν = measured_frequency(k, u)
		Jinv = emissivity_invariant(obj, x4, ν)
		Ainv = absorption_invariant(obj, x4, ν)
		Δλ = Δz / kz
		Iinv = _step_invariant(Iinv, Jinv, Ainv, Δλ)
	end

	νcam^3 * Iinv
end

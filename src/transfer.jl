struct AccValue{WHAT, T}
	value::T
end
AccValue{WHAT}(value) where {WHAT} = AccValue{WHAT, typeof(value)}(value)

_init_acc(::Type{T}, ν0) where {T<:Tuple} = map(a -> _init_acc(a, ν0), fieldtypes(T))
_init_acc(::Type{Intensity}, ν0) = AccValue{Intensity}(0)
_init_acc(::Type{IntensityIQU}, ν0) = begin
	AccValue{IntensityIQU}(StokesIQU(0, 0, 0))
end
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
	# ForwardDiff trick for ∂/∂ν (then convert to ∂/∂lnν downstream):
	# represent a dimensionless scale s = 1 + ε/ν₀ so that d(ν·s)/dε|₀ = 1.
	s = ForwardDiff.Dual{DT}(one(ν0f), inv(ν0f))
	Iinv = ForwardDiff.Dual{DT}(zero(ν0f), zero(ν0f))
	AccSpectralIndex(Iinv, s, DT)
end


_postprocess_acc(acc::Tuple, ν, what::Tuple) = map((a, w) -> _postprocess_acc(a, ν, w), acc, what)
_postprocess_acc(Iinv::AccValue{Intensity}, ν, what::Intensity) = Iinv.value * ν^3
_postprocess_acc(Iinv::AccValue{IntensityIQU}, ν, what::IntensityIQU) = Iinv.value * ν^3
# Optical depth is dimensionless and (with 𝓐 ≡ α_ν·ν and τ = ∫𝓐 dλ) is already
# the quantity used in transfer; unlike intensity, it does not require a ν³ conversion.
_postprocess_acc(τ::AccValue{OpticalDepth}, ν, what::OpticalDepth) = τ.value
_postprocess_acc(acc::AccSpectralIndex, ν, what::SpectralIndex) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	# Frequency conversion (invariant → ordinary): I_ν = ν³ 𝓘.
	# Therefore dI_ν/dν = ν³ d𝓘/dν + 3ν² 𝓘.
	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	# Spectral index (radio convention): α(ν) ≡ d ln I_ν / d ln ν = (ν/I_ν)·(dI_ν/dν).
	return (ν / I0) * dI
end
_postprocess_acc(acc::AccSpectralIndex, ν, what::Tuple{Intensity,SpectralIndex}) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	# Same math as above, but also return I_ν itself.
	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	return I0, (ν / I0) * dI
end

const Δτ_THRESHOLD_LINEAR = 1e-2


@inline _integrate_ray_step(acc::AccValue{Intensity}, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc.value
	# Invariant transfer: 𝓘 ≡ I_ν/ν³, 𝓙 ≡ j_ν/ν², 𝓐 ≡ α_ν·ν
	# d𝓘/dλ = 𝓙 − 𝓐 𝓘, so over a step Δλ with constant coeffs:
	# Δτ = 𝓐Δλ and 𝓘_out = 𝓘_in e^(−Δτ) + (𝓙/𝓐)(1−e^(−Δτ)).
	Δτ = Ainv * Δλ
	@assert Δτ ≥ 0
	Iinv = if Δτ < Δτ_THRESHOLD_LINEAR
        Iinv + (Jinv - Ainv * Iinv) * Δλ
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return AccValue{Intensity}(Iinv)
end

@inline _rt_step_scalar(y, Jinv, Ainv, Δλ) = begin
	Δτ = Ainv * Δλ
	@assert Δτ ≥ 0
	if Δτ < Δτ_THRESHOLD_LINEAR
		return y + (Jinv - Ainv * y) * Δλ
	end
	E = exp(-Δτ)
	return y * E + (Jinv / Ainv) * (1 - E)
end

@inline _integrate_ray_step(acc::AccValue{IntensityIQU}, obj::AbstractSynchrotronMedium, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv_m, Ainv_m) = emissivity_absorption_polarized_invariant(obj, x4, k′)

	# Build the comoving camera screen basis and the field-aligned +Q axis.
	(n′, e1′, e2′) = comoving_screen_basis(u, k)
	B′ = magnetic_field(obj, x4)
	(_, e_perp) = linear_polarization_basis_from_B(n′, B′)
	R_cf = stokes_QU_rotation(e1′, e2′, e_perp)  # (Q,U)_field = R_cf * (Q,U)_cam
	R_fc = R_cf'  # (Q,U)_cam = R_fc * (Q,U)_field

	# Camera → field Stokes.
	S_cam = acc.value
	S_f = rotate_QU(R_cf, S_cam)

	# Field Stokes → normal-mode state y = (I_perp, I_par, U).
	y_old = modes_from_IQU(S_f)

	# add U so that we can handle everything consistently
	Jinv_plus = ModePerpParU(Jinv_m.perp, Jinv_m.par, 0)
	Ainv_plus = ModePerpParU(Ainv_m.perp, Ainv_m.par, (Ainv_m.perp + Ainv_m.par)/2)
	
	# Diagonal normal-mode update (constant coeffs over the step).
	y_new = _rt_step_scalar.(y_old, Jinv_plus, Ainv_plus, Δλ)

	# Normal-mode → field Stokes.
	S_f_out = stokes_IQU(y_new)

	# Field → camera Stokes.
	QU_cam_out = rotate_QU(R_fc, @swiz S_f_out.QU)
	S_cam_out = StokesIQU(S_f_out.I, QU_cam_out...)
	return AccValue{IntensityIQU}(S_cam_out)
end

@inline _integrate_ray_step(acc::AccValue{OpticalDepth}, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	Ainv = absorption_invariant(obj, x4, k′)

	# Optical depth accumulation uses the invariant absorption 𝓐 = α_ν·ν:
	# τ = ∫ 𝓐 dλ, discretized as Σ (𝓐 Δλ).
	# The frame/Doppler factor is already inside 𝓐Δλ via ν′ = −k⋅u and Δλ = Δz/kᶻ:
	# for k = ν_obs(1,0,0,1), Δτ = (α′_ν ν′)(Δz/ν_obs) = α′_ν (Δz/δ), δ ≡ ν_obs/ν′.
	Δτ = Ainv * Δλ
	return AccValue{OpticalDepth}(acc.value + Δτ)
end

@inline _integrate_ray_step(acc::AccSpectralIndex, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	# Spectral index via AD: scale ν by s and rescale affine step.
	# Since k ∝ ν, scaling ν→ν·s implies λ-steps scale as Δλ→Δλ/s.
	k′ = lorentz_unboost(u, k * acc.s)
	Δλ′ = Δλ / acc.s
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc.Iinv
	Δτ = Ainv * Δλ′
	Δτ0 = ForwardDiff.value(Δτ)
	@assert Δτ0 ≥ 0
	Iinv = if Δτ0 < Δτ_THRESHOLD_LINEAR
		Iinv + (Jinv - Ainv * Iinv) * Δλ′
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return AccSpectralIndex(Iinv, acc.s, acc.DT)
end

_integrate_ray_step(acc::Tuple{AccValue{Intensity}, AccValue{OpticalDepth}}, obj, x4, k, Δλ) = begin
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc[1].value
	Δτ = Ainv * Δλ
	@assert Δτ ≥ 0
	Iinv = if Δτ < Δτ_THRESHOLD_LINEAR
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
	# Ray parameterization: k ∥ (1,0,0,1) so events satisfy x(z) = x₀ + z·(1,0,0,1).
	# With k = ν_cam·(1,0,0,1), we have kᶻ = ν_cam and Δλ = Δz / kᶻ.
    k1 = k / kz

	acc = if isempty(seg)
		z = leftendpoint(seg)
		_integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + z * k1, k, zero(float(z)))
	else
		# zs = range(seg, ray.nz)  # using StepRangeLen constructor directly is faster
		zs = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), ray.nz)
		Δz = step(zs)
		Δλ = Δz / kz
		acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + first(zs) * k1, k, Δλ)
		for z in zs[2:end]
			# Spacetime stepping along the null ray: x(z) = x₀ + z·k̂, with k̂ = (1,0,0,1).
			x = ray.x0 + z * k1
			acc = _integrate_ray_step(acc, obj, x, k, Δλ)
		end
		acc
	end

	ν = frequency(ray)
	return _postprocess_acc(acc, ν, what)
end


"""
	ray_contribution_profile(obj, ray) -> NamedTuple

Compute the *per-ray* contribution of each ray step (indexed by the internal z-grid)
to the final observed intensity at the camera.

This is a diagnostic helper for questions like:

- "How much does a point at depth z contribute to the current image pixel?"

Mathematically, for each step i the returned `dIinv_to_obs[i]` is the invariant source
term for that step, attenuated by the optical depth *in front of it* (toward the observer):

    dIinv_to_obs[i] = dIinv_source[i] * exp(-τ_front[i])

The discretization intentionally matches `integrate_ray` (same z-grid, same Δλ, and the
same linear-vs-exact handling for small Δτ).

Returns a NamedTuple with:

- `zs`: z-grid used for stepping
- `Δτ`: per-step invariant optical-depth increments (same indexing as `zs`)
- `τ_front`: cumulative optical depth in front of each step (τ from next step to exit)
- `dIinv_source`: per-step invariant source contribution before front attenuation
- `dIinv_to_obs`: per-step invariant contribution to the final pixel
- `dIν_to_obs`: same as `dIinv_to_obs`, converted to ordinary intensity via ν³
"""
ray_contribution_profile(obj::AbstractMedium, ray::RayZ) = begin
	seg = z_interval(obj, ray)
	k = ray.k
	kz = k.z
	@assert k == SVector(kz, 0, 0, kz)
	k1 = k / kz

	νobs = frequency(ray)
	FT = float(eltype(ray.x0))

	# Same z-grid and step sizing as integrate_ray.
	zs = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), isempty(seg) ? 0 : ray.nz)
	Δz = step(zs)
	Δλ = Δz / kz

	profile₁ = map(StructArray(; z=zs)) do (;z)
		x4 = ray.x0 + z * k1
		u = four_velocity(obj, x4)
		k′ = lorentz_unboost(u, k)
		(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)
		Δτ = Ainv * Δλ
		@assert Δτ ≥ 0
		dIinv_source = if Δτ < Δτ_THRESHOLD_LINEAR
			Jinv * Δλ
		else
			( Jinv / Ainv ) * (1 - exp(-Δτ))
		end
		(; z, x4, Δτ, Δz, dIinv_source)
	end

	# τ_front[i] = Σ_{j>i} Δτ[j]
	profile₂ = @insert profile₁.τ_front = reverse(cumsum(reverse(profile₁.Δτ))) .- profile₁.Δτ

	profile₃ = @insert profile₂.dIν_to_obs =
		map(profile₂) do (; dIinv_source, τ_front)
			dIinv_source * exp(-τ_front) * νobs^3
		end

	return profile₃
end

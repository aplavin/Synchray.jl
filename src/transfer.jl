struct AccValue{WHAT, T}
	value::T
end
AccValue{WHAT}(value) where {WHAT} = AccValue{WHAT, typeof(value)}(value)

_init_acc(::Type{T}, ν0) where {T<:Tuple} = map(a -> _init_acc(a, ν0), fieldtypes(T))
_init_acc(::Type{Intensity}, ν0) = AccValue{Intensity}(0)
_init_acc(::Type{IntensityIQU}, ν0) = AccValue{IntensityIQU}(StokesIQU(0, 0, 0))
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

const Δτ_THRESHOLD_LINEAR = 0.01f0


@inline _integrate_ray_step(acc::AccValue{Intensity}, obj, x4, ray, Δλ) = begin
	k = ray.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc.value
	# Invariant transfer: 𝓘 ≡ I_ν/ν³, 𝓙 ≡ j_ν/ν², 𝓐 ≡ α_ν·ν
	# d𝓘/dλ = 𝓙 − 𝓐 𝓘, so over a step Δλ with constant coeffs:
	# Δτ = 𝓐Δλ and 𝓘_out = 𝓘_in e^(−Δτ) + (𝓙/𝓐)(1−e^(−Δτ)).
	Δτ = Ainv * Δλ
	@boundscheck @assert Δτ ≥ 0
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
	@boundscheck @assert Δτ ≥ 0
	if Δτ < Δτ_THRESHOLD_LINEAR
		return y + (Jinv - Ainv * y) * Δλ
	end
	E = exp(-Δτ)
	return y * E + (Jinv / Ainv) * (1 - E)
end

@inline _integrate_ray_step(acc::AccValue{IntensityIQU}, obj::AbstractMedium, x4, ray, Δλ) = begin
	k = ray.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv_m, Ainv_m, B′) = emissivity_absorption_polarized_invariant(obj, x4, k′)

	# Build the comoving camera screen basis and the field-aligned +Q axis.
	(n′, e1′, e2′) = comoving_screen_basis(u, ray)
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

@inline _integrate_ray_step(acc::AccValue{OpticalDepth}, obj, x4, ray, Δλ) = begin
	k = ray.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	Ainv = absorption_invariant(obj, x4, k′)

	# Optical depth accumulation uses the invariant absorption 𝓐 = α_ν·ν:
	# τ = ∫ 𝓐 dλ, discretized as Σ (𝓐 Δλ).
	Δτ = Ainv * Δλ
	return AccValue{OpticalDepth}(acc.value + Δτ)
end

@inline _integrate_ray_step(acc::AccSpectralIndex, obj, x4, ray, Δλ) = begin
	k = ray.k
	u = four_velocity(obj, x4)
	# Spectral index via AD: scale ν by s and rescale affine step.
	# Since k ∝ ν, scaling ν→ν·s implies λ-steps scale as Δλ→Δλ/s.
	k′ = lorentz_unboost(u, k * acc.s)
	Δλ′ = Δλ / acc.s
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc.Iinv
	Δτ = Ainv * Δλ′
	Δτ0 = ForwardDiff.value(Δτ)
	@boundscheck @assert Δτ0 ≥ 0
	Iinv = if Δτ0 < Δτ_THRESHOLD_LINEAR
		Iinv + (Jinv - Ainv * Iinv) * Δλ′
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return AccSpectralIndex(Iinv, acc.s, acc.DT)
end

_integrate_ray_step(acc::Tuple{AccValue{Intensity}, AccValue{OpticalDepth}}, obj, x4, ray, Δλ) = begin
	k = ray.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)

	Iinv = acc[1].value
	Δτ = Ainv * Δλ
	@boundscheck @assert Δτ ≥ 0
	Iinv = if Δτ < Δτ_THRESHOLD_LINEAR
        Iinv + (Jinv - Ainv * Iinv) * Δλ
	else
		E = exp(-Δτ)
		Iinv * E + (Jinv / Ainv) * (1 - E)
	end
	return (AccValue{Intensity}(Iinv),
			AccValue{OpticalDepth}(acc[2].value + Δτ))
end


integrate_ray(obj::AbstractMedium, ray::Ray, what=Intensity()) = begin
	seg = z_interval(obj, ray)

	k = ray.k
	ν = frequency(ray)
	k1 = direction4(k)

	acc = if isempty(seg)
		s = leftendpoint(seg)
		_integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + s * k1, ray, zero(float(s)))
	else
		# zs = range(seg, ray.nz)  # using StepRangeLen constructor directly is faster
		@boundscheck @assert ray.nz ≥ 0  # to convince compiler that StepRangeLen won't error, necessary for running in Metal
		ss = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), ray.nz)
		Δs = step(ss)
		Δλ = Δs / ν
		acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + first(ss) * k1, ray, Δλ)
		for s in Iterators.drop(ss, 1)  # cannot use ss[2:end] because StepRangeLen indexing can throw, failing in Metal
			x = ray.x0 + s * k1
			acc = _integrate_ray_step(acc, obj, x, ray, Δλ)
		end
		acc
	end

	return _postprocess_acc(acc, ν, what)
end


# --- CombinedMedium integration ---
@inline _integrate_segment(acc, obj, ray, seg) = begin
	isempty(seg) && return acc
	k = ray.k
	ν = frequency(ray)
	k1 = direction4(k)
	@boundscheck @assert ray.nz ≥ 0
	ss = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), ray.nz)
	Δs = step(ss)
	Δλ = Δs / ν
	for s in ss
		x = ray.x0 + s * k1
		acc = _integrate_ray_step(acc, obj, x, ray, Δλ)
	end
	return acc
end

# N=1: delegate directly, zero overhead
integrate_ray(cm::CombinedMedium{<:Tuple{Any}}, ray::Ray, what=Intensity()) = integrate_ray(cm.objects[1], ray, what)

# General: integrate objects in tuple order (must be back-to-front)
integrate_ray(cm::CombinedMedium, ray::Ray, what=Intensity()) = begin
	k = ray.k
	ν = frequency(ray)
	k1 = direction4(k)

	# Type-promoting zero-step (same pattern as existing integrate_ray)
	obj1 = first(cm.objects)
	seg1 = z_interval(obj1, ray)
	s_init = leftendpoint(seg1)
	acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj1, ray.x0 + s_init * k1, ray, zero(float(s_init)))

	acc = _integrate_combined_recursive(acc, cm.objects, ray)

	return _postprocess_acc(acc, ν, what)
end

@inline _integrate_combined_recursive(acc, ::Tuple{}, ray) = acc
@inline _integrate_combined_recursive(acc, objs::Tuple, ray) = begin
	obj = first(objs)
	seg = z_interval(obj, ray)
	acc = _integrate_segment(acc, obj, ray, seg)
	_integrate_combined_recursive(acc, Base.tail(objs), ray)
end


"""
	ray_contribution_profile(obj, ray) -> NamedTuple

Compute the *per-ray* contribution of each ray step to the final observed intensity.

Mathematically, for each step i the returned `dIinv_to_obs[i]` is the invariant source
term for that step, attenuated by the optical depth *in front of it* (toward the observer):

    dIinv_to_obs[i] = dIinv_source[i] * exp(-τ_front[i])

Returns a StructArray with columns:
- `zs`: z-grid used for stepping
- `Δτ`: per-step invariant optical-depth increments (same indexing as `zs`)
- `τ_front`: cumulative optical depth in front of each step (τ from next step to exit)
- `dIinv_source`: per-step invariant source contribution before front attenuation
- `dIinv_to_obs`: per-step invariant contribution to the final pixel
- `dIν_to_obs`: same as `dIinv_to_obs`, converted to ordinary intensity via ν³
"""
ray_contribution_profile(obj::AbstractMedium, ray::Ray) = begin
	seg = z_interval(obj, ray)
	k = ray.k
	ν = frequency(ray)
	k1 = direction4(k)

	νobs = frequency(ray)

	# Same grid and step sizing as integrate_ray.
	ss = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), isempty(seg) ? 0 : ray.nz)
	Δs = step(ss)
	Δλ = Δs / ν

	profile₁ = map(StructArray(; s=ss)) do (;s)
		x4 = ray.x0 + s * k1
		u = four_velocity(obj, x4)
		k′ = lorentz_unboost(u, k)
		(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′)
		Δτ = Ainv * Δλ
		@boundscheck @assert Δτ ≥ 0
		dIinv_source = if Δτ < Δτ_THRESHOLD_LINEAR
			Jinv * Δλ
		else
			( Jinv / Ainv ) * (1 - exp(-Δτ))
		end
		(; s, x4, Δτ, Δs, dIinv_source)
	end

	# τ_front[i] = Σ_{j>i} Δτ[j]
	profile₂ = @insert profile₁.τ_front = reverse(cumsum(reverse(profile₁.Δτ))) .- profile₁.Δτ

	profile₃ = @insert profile₂.dIν_to_obs =
		map(profile₂) do (; dIinv_source, τ_front)
			dIinv_source * exp(-τ_front) * νobs^3
		end

	return profile₃
end


"""
	ray_contribution_profile_IQU(obj, ray) -> StructArray

Compute the *per-ray* contribution of each ray step to the final observed Stokes vector
(I, Q, U) at the camera, accounting for polarization-dependent emission, absorption, and
basis rotations.

This extends `ray_contribution_profile` to full Stokes polarization.

Mathematical equivalence:
    sum(profile.dIν_to_obs) ≈ render(ray, obj, IntensityIQU())

Returns a StructArray with the same structure as `ray_contribution_profile`, but with
`dIν_to_obs` as a StokesIQU vector instead of scalar.
"""
ray_contribution_profile_IQU(obj::AbstractMedium, ray::Ray) = begin
	seg = z_interval(obj, ray)
	k = ray.k
	ν = frequency(ray)
	k1 = direction4(k)

	νobs = frequency(ray)

	# Same grid and step sizing as integrate_ray
	ss = StepRangeLen(leftendpoint(seg), width(seg) / (ray.nz - 1), isempty(seg) ? 0 : ray.nz)
	Δs = step(ss)
	Δλ = Δs / ν

	# Step 1: Forward pass - collect geometry, coefficients, and basis rotations
	profile₁ = map(StructArray(; s=ss)) do (;s)
		x4 = ray.x0 + s * k1
		u = four_velocity(obj, x4)
		k′ = lorentz_unboost(u, k)
		ν′ = frequency(k′)

		# Get polarized coefficients
		(j_modes, α_modes, B′) = emissivity_absorption_polarized(obj, x4, k′)
		Jinv_modes = j_modes / (ν′^2)
		Ainv_modes = α_modes * ν′

		# Optical depth per step (mode-wise)
		Δτ_modes = Ainv_modes * Δλ

		# Get basis rotation: camera → field
		(n′, e1′, e2′) = comoving_screen_basis(u, ray)
		(e_par, e_perp) = linear_polarization_basis_from_B(n′, B′)
		R_camera_to_field = stokes_QU_rotation(e1′, e2′, e_perp)

		(; s, x4, Δs, Jinv_modes, Ainv_modes, Δτ_modes, R_camera_to_field)
	end

	# Early return for empty profile (ray misses geometry)
	if isempty(profile₁)
		return StructArray((; s=eltype(ss)[], x4=eltype(profile₁.x4)[], Δs=eltype(profile₁.Δs)[], dIν_to_obs=StokesIQU{Float64}[]))
	end

	# Step 2: Cumulative front optical depth (mode-wise)
	Δτ_perp_vec = [step.Δτ_modes.perp for step in profile₁]
	Δτ_par_vec = [step.Δτ_modes.par for step in profile₁]
	τ_front_perp = reverse(cumsum(reverse(Δτ_perp_vec))) .- Δτ_perp_vec
	τ_front_par = reverse(cumsum(reverse(Δτ_par_vec))) .- Δτ_par_vec

	# Step 3-5: Emission, attenuation, rotation to camera frame
	dIν_to_obs = map(eachindex(profile₁)) do i
		step = profile₁[i]

		# Emission term accounting for step's internal opacity
		dSinv_emit_modes = if all(step.Δτ_modes .< Δτ_THRESHOLD_LINEAR)
			step.Jinv_modes * Δλ
		else
			ModePerpPar(
				step.Δτ_modes.perp < Δτ_THRESHOLD_LINEAR ?
					step.Jinv_modes.perp * Δλ :
					(step.Jinv_modes.perp / step.Ainv_modes.perp) * (1 - exp(-step.Δτ_modes.perp)),
				step.Δτ_modes.par < Δτ_THRESHOLD_LINEAR ?
					step.Jinv_modes.par * Δλ :
					(step.Jinv_modes.par / step.Ainv_modes.par) * (1 - exp(-step.Δτ_modes.par))
			)
		end

		# Convert modes to field-frame Stokes
		I_emit = dSinv_emit_modes.perp + dSinv_emit_modes.par
		Q_emit = dSinv_emit_modes.perp - dSinv_emit_modes.par
		U_emit = zero(I_emit)

		# Front material attenuation (mode-wise)
		E_perp = exp(-τ_front_perp[i])
		E_par = exp(-τ_front_par[i])

		I_perp_emit = (I_emit + Q_emit) / 2
		I_par_emit = (I_emit - Q_emit) / 2

		I_perp_atten = I_perp_emit * E_perp
		I_par_atten = I_par_emit * E_par

		τ_U_front = (τ_front_perp[i] + τ_front_par[i]) / 2
		U_atten = U_emit * exp(-τ_U_front)

		I_atten = I_perp_atten + I_par_atten
		Q_atten = I_perp_atten - I_par_atten

		# Rotate to camera frame
		R_field_to_camera = step.R_camera_to_field'
		QU_camera = R_field_to_camera * SVector(Q_atten, U_atten)

		# Convert to observer frequency
		StokesIQU(I_atten, QU_camera[1], QU_camera[2]) * νobs^3
	end

	return StructArray((; profile₁.s, profile₁.x4, profile₁.Δs, dIν_to_obs))
end

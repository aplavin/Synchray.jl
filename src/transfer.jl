struct AccValue{WHAT, T}
	value::T
end
AccValue{WHAT}(value) where {WHAT} = AccValue{WHAT, typeof(value)}(value)

_init_acc(::Type{T}, ν0) where {T<:Tuple} = map(a -> _init_acc(a, ν0), fieldtypes(T))
# Float (not Int) init: Int seed makes `acc` type-unstable → Metal-throwing union-split
_init_acc(::Type{Intensity}, ν0) = AccValue{Intensity}(zero(float(ν0)))
_init_acc(::Type{IntensityIQU}, ν0) = let z = zero(float(ν0)); AccValue{IntensityIQU}(StokesIQU(z, z, z)) end
_init_acc(::Type{OpticalDepth}, ν0) = AccValue{OpticalDepth}(zero(float(ν0)))
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
# The accumulator carries 𝓘̂ = I_ν·(ν_ref/ν)³, which equals the physical I_ν at the observer.
_postprocess_acc(Iinv::AccValue{Intensity}, ν, what::Intensity) = Iinv.value
_postprocess_acc(Iinv::AccValue{IntensityIQU}, ν, what::IntensityIQU) = Iinv.value
# Optical depth is dimensionless and (with 𝓐 ≡ α_ν·ν and τ = ∫𝓐 dλ) is already
# the quantity used in transfer; unlike intensity, it does not require a ν³ conversion.
_postprocess_acc(τ::AccValue{OpticalDepth}, ν, what::OpticalDepth) = τ.value
_postprocess_acc(acc::AccSpectralIndex, ν, what::SpectralIndex) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	# Radio-convention spectral index α(ν) ≡ d ln I_ν / d ln ν = (ν/I_ν)·(dI_ν/dν).
	# I0 and dI are I_ν and dI_ν/dν up to a common constant factor, which cancels in the ratio.
	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	return (ν / I0) * dI
end
_postprocess_acc(acc::AccSpectralIndex, ν, what::Tuple{Intensity,SpectralIndex}) = begin
	(;Iinv, DT) = acc
	Iinv0 = ForwardDiff.value(DT, Iinv)
	dIinv = ForwardDiff.extract_derivative(DT, Iinv)

	# Spectral index as above, plus the intensity itself: Iinv0 already equals physical I_ν.
	I0 = Iinv0 * ν^3
	dI = dIinv * ν^3 + Iinv0 * (3 * ν^2)
	return Iinv0, (ν / I0) * dI
end

const Δτ_THRESHOLD_LINEAR = 0.01f0


# One transfer step at x4 over Δλ. A `CombinedMedium` applies the step of each component containing
# x4, in object order (overlaps: sequential splitting, exact as Δλ → 0). Per-component steps keep
# `acc` type-stable for heterogeneous components. `@propagate_inbounds`: a caller's `@inbounds`
# must reach the `@boundscheck`s of the single-medium steps.
Base.@propagate_inbounds _integrate_ray_step(acc, obj, x4, st::RayState, Δλ) = _integrate_ray_step_single(acc, obj, x4, st, Δλ)
Base.@propagate_inbounds _integrate_ray_step(acc, cm::CombinedMedium, x4, st::RayState, Δλ) = _step_each_inside(acc, cm.objects, x4, st, Δλ)

# Heterogeneous tuple: recursive unrolling (no dynamic dispatch, GPU-safe).
@inline _step_each_inside(acc, ::Tuple{}, x4, st, Δλ) = acc
Base.@propagate_inbounds _step_each_inside(acc, objs::Tuple, x4, st, Δλ) = begin
	obj = first(objs)
	acc = is_inside(obj, x4) ? _integrate_ray_step(acc, obj, x4, st, Δλ) : acc
	_step_each_inside(acc, Base.tail(objs), x4, st, Δλ)
end
# Homogeneous tuples and vectors: a loop, so the step body is compiled once, not per component.
Base.@propagate_inbounds _step_each_inside(acc, objs::NTuple{N,T} where {N,T}, x4, st, Δλ) = _step_each_inside_loop(acc, objs, x4, st, Δλ)
Base.@propagate_inbounds _step_each_inside(acc, objs::AbstractVector, x4, st, Δλ) = _step_each_inside_loop(acc, objs, x4, st, Δλ)
Base.@propagate_inbounds _step_each_inside_loop(acc, objs, x4, st, Δλ) = begin
	for obj in objs
		is_inside(obj, x4) && (acc = _integrate_ray_step(acc, obj, x4, st, Δλ))
	end
	acc
end


@inline _integrate_ray_step_single(acc::AccValue{Intensity}, obj, x4, st, Δλ) = begin
	k = st.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′, st.ν_ref)

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
	if Δτ < Δτ_THRESHOLD_LINEAR
		return y + (Jinv - Ainv * y) * Δλ
	end
	E = exp(-Δτ)
	return y * E + (Jinv / Ainv) * (1 - E)
end

@inline _integrate_ray_step_single(acc::AccValue{IntensityIQU}, obj::AbstractMedium, x4, st, Δλ) = begin
	k = st.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv_m, Ainv_m, B′) = emissivity_absorption_polarized_invariant(obj, x4, k′, st.ν_ref)

	# Build the comoving camera screen basis and the field-aligned +Q axis.
	(n′, e1′, e2′) = comoving_screen_basis(u, st)
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

@inline _integrate_ray_step_single(acc::AccValue{OpticalDepth}, obj, x4, st, Δλ) = begin
	k = st.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	Ainv = absorption_invariant(obj, x4, k′)

	# Optical depth accumulation uses the invariant absorption 𝓐 = α_ν·ν:
	# τ = ∫ 𝓐 dλ, discretized as Σ (𝓐 Δλ).
	Δτ = Ainv * Δλ
	return AccValue{OpticalDepth}(acc.value + Δτ)
end

@inline _integrate_ray_step_single(acc::AccSpectralIndex, obj, x4, st, Δλ) = begin
	k = st.k
	u = four_velocity(obj, x4)
	# Spectral index via AD: scale ν by s and rescale affine step.
	# Since k ∝ ν, scaling ν→ν·s implies λ-steps scale as Δλ→Δλ/s.
	k′ = lorentz_unboost(u, k * acc.s)
	Δλ′ = Δλ / acc.s
	# st.ν_ref is the unperturbed reference frequency, not the AD-scaled k*acc.s.
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′, st.ν_ref)

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

_integrate_ray_step_single(acc::Tuple{AccValue{Intensity}, AccValue{OpticalDepth}}, obj, x4, st, Δλ) = begin
	k = st.k
	u = four_velocity(obj, x4)
	k′ = lorentz_unboost(u, k)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, k′, st.ν_ref)

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


# Each midpoint represents one full cell; the nz cell widths sum to the clipped
# segment length. Shared by rendering and contribution profiles.
@inline _segment_midpoints(seg, nz) = begin
	@boundscheck @assert nz ≥ 1 "nz must be ≥ 1"
	Δs = width(seg) / nz
	StepRangeLen(leftendpoint(seg) + Δs / 2, Δs, isempty(seg) ? 0 : nz)
end


integrate_ray(obj::AbstractMedium, ray::Ray, what=Intensity()) = begin
	seg = z_interval_clipped(obj, ray)

	ν = frequency(ray)
	k1 = direction4(ray)
	st = RayState(ray)

	acc = if isempty(seg)
		s = leftendpoint(seg)
		_integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + s * k1, st, zero(float(s)))
	else
		ss = _segment_midpoints(seg, ray.nz)
		Δs = step(ss)
		Δλ = Δs / ν
		acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + first(ss) * k1, st, Δλ)
		for s in Iterators.drop(ss, 1)  # cannot use ss[2:end] because StepRangeLen indexing can throw, failing in Metal
			x = ray.x0 + s * k1
			acc = _integrate_ray_step(acc, obj, x, st, Δλ)
		end
		acc
	end

	return _postprocess_acc(acc, ν, what)
end


# --- CombinedMedium integration ---
# XXX: @inline hurts CPU performance a bit with 3+ objects, but without it GPU is slow
@inline _integrate_segment(acc, obj, ray, seg) = begin
	isempty(seg) && return acc
	ν = frequency(ray)
	k1 = direction4(ray)
	st = RayState(ray)
	ss = _segment_midpoints(seg, ray.nz)
	Δs = step(ss)
	Δλ = Δs / ν
	for s in ss
		x = ray.x0 + s * k1
		acc = _integrate_ray_step(acc, obj, x, st, Δλ)
	end
	return acc
end

# N=1: delegate directly, zero overhead
integrate_ray(cm::CombinedMedium{<:Tuple{Any}}, ray::Ray, what=Intensity()) = integrate_ray(cm.objects[1], ray, what)

# Tuple: recursive unrolling for type stability and GPU
integrate_ray(cm::CombinedMedium{<:Tuple}, ray::Ray, what=Intensity()) = begin
	ν = frequency(ray)
	k1 = direction4(ray)

	# Type-promoting zero-step (same pattern as existing integrate_ray)
	obj1 = first(cm.objects)
	seg1 = z_interval_clipped(obj1, ray)
	s_init = leftendpoint(seg1)
	acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj1, ray.x0 + s_init * k1, RayState(ray), zero(float(s_init)))

	acc = _integrate_combined_recursive(acc, cm.objects, ray)

	return _postprocess_acc(acc, ν, what)
end

@inline _integrate_combined_recursive(acc, ::Tuple{}, ray) = acc
@inline _integrate_combined_recursive(acc, objs::Tuple, ray) = begin
	obj = first(objs)
	seg = z_interval_clipped(obj, ray)
	acc = _integrate_segment(acc, obj, ray, seg)
	_integrate_combined_recursive(acc, Base.tail(objs), ray)
end

# AbstractVector: loop-based
integrate_ray(cm::CombinedMedium{<:AbstractVector}, ray::Ray, what=Intensity()) = begin
	ν = frequency(ray)
	k1 = direction4(ray)

	obj1 = first(cm.objects)
	seg1 = z_interval_clipped(obj1, ray)
	s_init = leftendpoint(seg1)
	acc = _integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj1, ray.x0 + s_init * k1, RayState(ray), zero(float(s_init)))

	for obj in cm.objects
		seg = z_interval_clipped(obj, ray)
		acc = _integrate_segment(acc, obj, ray, seg)
	end

	return _postprocess_acc(acc, ν, what)
end


# Source of one constant-coefficient step as seen at its exit: ∫₀^Δλ 𝓙̂·e^(−𝓐·(Δλ−λ)) dλ.
@inline _step_source(Jinv, Ainv, Δλ) = begin
	Δτ = Ainv * Δλ
	Δτ < Δτ_THRESHOLD_LINEAR ? Jinv * Δλ : (Jinv / Ainv) * (1 - exp(-Δτ))
end

"""
	step_contribution(obj, x4, st::RayState, Δλ) -> (; Δτ, dIinv_source)

Optical depth `Δτ` of one transfer step of `obj` at `x4` over `Δλ`, and the step's source
`dIinv_source` (physical intensity scale) before attenuation by the medium in front of it.
Uses the pointwise coefficients of `obj`.
"""
@inline step_contribution(obj, x4, st::RayState, Δλ) = begin
	u = four_velocity(obj, x4)
	(Jinv, Ainv) = emissivity_absorption_invariant(obj, x4, lorentz_unboost(u, st.k), st.ν_ref)
	Δτ = Ainv * Δλ
	@boundscheck @assert Δτ ≥ 0
	(; Δτ, dIinv_source = _step_source(Jinv, Ainv, Δλ))
end

"""
	attenuate_to_observer(profile) -> StructArray

Add the columns `τ_front` (optical depth in front of each step, toward the observer) and
`dIν_to_obs = dIinv_source·exp(−τ_front)` (the step's contribution to the observed intensity)
to a per-step profile with columns `Δτ` and `dIinv_source`, ordered from source to observer.
"""
attenuate_to_observer(profile) = begin
	profile₁ = @insert profile.τ_front = _front_optical_depth(profile.Δτ)
	@insert profile₁.dIν_to_obs = map(profile₁) do (; dIinv_source, τ_front)
		dIinv_source * exp(-τ_front)
	end
end

# τ_front[i] = Σ_{j>i} Δτ[j]
_front_optical_depth(Δτ) = reverse(cumsum(reverse(Δτ))) .- Δτ


"""
	ray_contribution_profile(obj, ray) -> StructArray

Compute the *per-ray* contribution of each ray step to the final observed intensity.

Mathematically, for each step i the returned `dIν_to_obs[i]` is the source term for that
step, attenuated by the optical depth *in front of it* (toward the observer):

    dIν_to_obs[i] = dIinv_source[i] * exp(-τ_front[i])

Returns a StructArray with columns:
- `s`: ray-parameter cell midpoints used for stepping
- `x4`: cell midpoint 4-positions
- `Δτ`: per-cell optical-depth increments (same indexing as `s`)
- `Δs`: cell widths
- `dIinv_source`: per-step source contribution (physical intensity scale) before front attenuation
- `τ_front`: cumulative optical depth in front of each step (τ from next step to exit)
- `dIν_to_obs`: per-step contribution to the final pixel (physical specific intensity)
"""
ray_contribution_profile(obj::AbstractMedium, ray::Ray) = begin
	seg = z_interval_clipped(obj, ray)
	k1 = direction4(ray)
	st = RayState(ray)

	# Same grid and step sizing as integrate_ray.
	ss = _segment_midpoints(seg, ray.nz)
	Δs = step(ss)
	Δλ = Δs / frequency(ray)

	profile = map(StructArray(; s=ss)) do (; s)
		x4 = ray.x0 + s * k1
		(; Δτ, dIinv_source) = step_contribution(obj, x4, st, Δλ)
		(; s, x4, Δτ, Δs, dIinv_source)
	end
	return attenuate_to_observer(profile)
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
	seg = z_interval_clipped(obj, ray)
	k1 = direction4(ray)
	st = RayState(ray)

	# Same grid and step sizing as integrate_ray.
	ss = _segment_midpoints(seg, ray.nz)
	Δs = step(ss)
	Δλ = Δs / frequency(ray)

	# Per step: mode-wise optical depth and source, and the camera → field (Q,U) rotation.
	profile = map(StructArray(; s=ss)) do (; s)
		x4 = ray.x0 + s * k1
		u = four_velocity(obj, x4)
		(Jinv_modes, Ainv_modes, B′) = emissivity_absorption_polarized_invariant(obj, x4, lorentz_unboost(u, st.k), st.ν_ref)
		(n′, e1′, e2′) = comoving_screen_basis(u, st)
		(_, e_perp) = linear_polarization_basis_from_B(n′, B′)
		(; s, x4, Δs,
		   Δτ_modes = Ainv_modes * Δλ,
		   dSinv_modes = _step_source.(Jinv_modes, Ainv_modes, Δλ),
		   R_camera_to_field = stokes_QU_rotation(e1′, e2′, e_perp))
	end

	# Each mode is attenuated by its own front optical depth; the emitted field-frame U is zero.
	dIν_to_obs = map(profile.dSinv_modes, _front_optical_depth(profile.Δτ_modes), profile.R_camera_to_field) do dS, τ_front, R
		(; I, Q) = stokes_IQ(dS .* exp.(-τ_front))
		StokesIQU(I, (R' * SVector(Q, zero(Q)))...)
	end
	return StructArray((; profile.s, profile.x4, profile.Δs, dIν_to_obs))
end

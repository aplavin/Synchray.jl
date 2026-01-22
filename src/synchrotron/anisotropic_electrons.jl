"""
		AnisotropicPowerLawElectrons(; p, η=1, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)

Power-law synchrotron electron model for Stokes-I with *gyrotropic* pitch-angle anisotropy.

This implements the simple prescription used by Tsunetoe et al. (2025) following
Dexter (2016) / Melrose (1971): for ordered fields, the emissivity and absorption
coefficients scale as

	j_I(θ_{Bn}) = φ(θ_{Bn}) · j_{I,iso}(θ_{Bn}),
	α_I(θ_{Bn}) = φ(θ_{Bn}) · α_{I,iso}(θ_{Bn}),

with

	φ(θ) = P(p,η)^{-1} · [1 + (η-1) cos²θ]^{-p/2},
	P(p,η) = ∫₀¹ dμ [1 + (η-1) μ²]^{-p/2}.

Here θ_{Bn} is the angle between the comoving magnetic-field direction and the
comoving photon direction.

Semantics of η:

- η = 1: isotropic (recovers `IsotropicPowerLawElectrons`)
- η < 1: electrons concentrated along the field (small pitch angles)
- η > 1: electrons concentrated perpendicular to the field

Notes / current scope:

- This model is only meaningful when the magnetic field has a direction (ordered vector
  or `TangledOrderedMixture`).
- For `FullyTangled(|B′|)`, anisotropy is rejected unless η==1.
"""
struct AnisotropicPowerLawElectrons{Tp,Tη,Tγ,TC,Tavg,TP}
	p::Tp
	η::Tη
	γmin::Tγ
	γmax::Tγ
	Cj_ordered::TC
	Ca_ordered::TC
	sinavg_j::Tavg
	sinavg_a::Tavg
	Pnorm::TP
end

@unstable prepare_for_computations(model::AnisotropicPowerLawElectrons) = @modify(FixedExponent, model.p)

# In the isotropic model, θ_{Bn} enters only through B_perp = |B'| sinθ_{Bn}.
# Here we additionally apply an anisotropy weight φ(θ_{Bn}), but since φ depends only on cos²θ,
# we take cos²θ directly.
@inline _phi_theta(model::AnisotropicPowerLawElectrons, cos2θ) = begin
	(;η, p, Pnorm) = model
	cos2 = clamp(cos2θ, 0, 1)
	return (1 + (η - 1) * cos2)^(-_half(p)) / Pnorm
end

# Constructor mirrors `IsotropicPowerLawElectrons`:
# - same ordered-field normalization for (Cj_ordered, Ca_ordered)
# - same sin-averages stored for compatibility with `TangledOrderedMixture`
# Difference: we also precompute P(p,η) and store it as `Pnorm` so coefficient evaluation stays cheap.
function AnisotropicPowerLawElectrons(; p, η=1, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)
	@assert p > 1 && γmin > 0 && γmax > γmin && η > 0

	qj = _half(p + 1)
	qa = _half(p + 2)
	sinavg_j = _avg_sin_pow(qj)
	sinavg_a = _avg_sin_pow(qa)

	if isnothing(Cj) || isnothing(Ca)
		@assert isnothing(Cj) && isnothing(Ca)
		K_per_ne = _K_per_ne(p, γmin, γmax)
		(c5, c6) = _synchrotron_c5_c6_ordered(p)
		Cj_ordered = c5 * K_per_ne
		Ca_ordered = c6 * K_per_ne
	else
		# Match `IsotropicPowerLawElectrons`: assume provided Cj/Ca were for tangled field.
		Cj_ordered = Cj / sinavg_j
		Ca_ordered = Ca / sinavg_a
	end

	Pnorm = let
		val, err = quadgk(μ -> (1 + (η - 1) * μ^2)^(-p / 2), 0, 1; rtol=1e-10)
		@assert err < 1e-5 * val
		val
	end
	return AnisotropicPowerLawElectrons(p, η, promote(γmin, γmax)..., promote(Cj_ordered, Ca_ordered)..., promote(sinavg_j, sinavg_a)..., Pnorm)
end

# Unlike the isotropic model, anisotropy is not defined for a directionless tangled field.
# We reject `FullyTangled` entirely rather than introducing an extra angle-averaged interpretation.
@inline _synchrotron_coeffs(model::AnisotropicPowerLawElectrons, n_e, field::FullyTangled, k′::FourFrequency) =
    error("AnisotropicPowerLawElectrons requires an ordered magnetic-field direction, FullyTangled is not supported")

# This corresponds to the ordered-vector `_synchrotron_coeffs` in `isotropic_electrons.jl`.
# Difference: multiply both coefficients by φ(θ_{Bn}; η), where η is stored on `model`
# (via `_phi_theta(model, cos2θ)` and cached `model.Pnorm`).
@inline _synchrotron_coeffs(model::AnisotropicPowerLawElectrons, n_e, b::SVector{3}, k′::FourFrequency) = let
	ν = frequency(k′)
	invν = inv(ν)
	(;p, Cj_ordered, Ca_ordered) = model

	n̂ = (@swiz k′.xyz) * invν
	b2 = dot(b, b)
	n2 = dot(n̂, n̂)
	dotbn = dot(b, n̂)
	Bperp = sqrt(max(b2 * n2 - dotbn^2, 0))
	cos2θ = dotbn^2 / (b2 * n2)
	φ = _phi_theta(model, cos2θ)

	# Previous version (kept for reference; arguably cleaner, but has another sqrt call):
	# B = norm(b)
	# Bperp = norm(cross(b, n̂))
	# cosθ = dot(b, n̂) / B
	# φ = _phi_theta(model, cosθ^2)

	B_over_ν = Bperp * invν
	common = B_over_ν^_half(p)
	j = Cj_ordered * n_e * common * sqrt(Bperp * ν) * φ
	α = Ca_ordered * n_e * common * Bperp * invν^2 * φ
	return j, α
end

# Unified TangledOrderedMixture method for both isotropic and anisotropic electrons (same pattern as polarization).
# This is defined here (in anisotropic_electrons.jl) because both types need to be loaded for the Union type.
@inline _synchrotron_coeffs(
	model::Union{IsotropicPowerLawElectrons, AnisotropicPowerLawElectrons},
	n_e,
	field::TangledOrderedMixture,
	k′::FourFrequency
) = let
	ν = frequency(k′)
	(;p, Cj_ordered, Ca_ordered, sinavg_j, sinavg_a) = model
	κ = field.kappa
	@assert κ ≥ 0

	b = field.b
	B = norm(b)
	invν = inv(ν)

	# Ordered viewing angle from the preferred direction.
	n = (@swiz k′.xyz) * invν
	@assert dot(n, n) ≈ 1
	sinθ = norm(cross(b, n)) / B
	sinθ = clamp(sinθ, 0, 1)
	cosθ = dot(b, n) / B

	# Minimal ordering model: mix between isotropic-direction average (κ=0) and fully ordered (κ→∞).
	f = κ == Inf ? one(float(κ)) : κ / (one(κ) + κ)

	qj = _half(p + StaticNum{1}())
	qa = _half(p + StaticNum{2}())

	# Dispatch to model-specific pitch-angle factor
	# (φ=1 for isotropic, φ(θ) for anisotropic)
	φ = _phi_theta(model, cosθ^2)

	# For anisotropic electrons: tangled component uses isotropic average,
	# ordered component includes anisotropy factor φ(θ_Bn).
	# For isotropic electrons: φ=1, so this reduces to the previous formula.
	Aj = muladd(f, φ * sinθ^qj - sinavg_j, sinavg_j)
	Aa = muladd(f, φ * sinθ^qa - sinavg_a, sinavg_a)

	B_over_ν = B * invν
	common = B_over_ν^_half(p)
	j = Cj_ordered * n_e * common * sqrt(B * ν) * Aj
	α = Ca_ordered * n_e * common * B * invν^2 * Aa
	return j, α
end

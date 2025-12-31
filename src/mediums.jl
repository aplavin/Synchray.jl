abstract type AbstractMedium end

abstract type AbstractSynchrotronMedium <: AbstractMedium end

# Primary API (what the raytracer should prefer calling): compute both comoving j_ν and α_ν
# (they typically share expensive intermediate computations).
#
# Radiative-transfer invariants used by the integrator:
#   𝓘 ≡ I_ν / ν^3
#   𝓙 ≡ j_ν / ν^2
#   𝓐 ≡ α_ν · ν
# so along an affine parameter λ:
#   d𝓘/dλ = 𝓙 - 𝓐 𝓘.
#
# This helper takes (j_ν, α_ν) at the frequency ν measured in the medium rest frame
# and returns the invariant pair (𝓙, 𝓐) used for integration.
@inline emissivity_absorption_invariant(obj::AbstractMedium, x4, ν) = begin
	(j, α) = emissivity_absorption(obj, x4, ν)
	return (j / (ν^2), α * ν)
end

# if a medium defines only `emissivity_absorption`, the generic `emissivity`/`absorption` wrappers below will work
# `absorption` is used in optical depth calculation, maybe can drop this in future?..
@inline emissivity_invariant(obj::AbstractMedium, x4, ν) = emissivity(obj, x4, ν) / (ν^2)
@inline absorption_invariant(obj::AbstractMedium, x4, ν) = absorption(obj, x4, ν) * ν
@inline emissivity(obj::AbstractMedium, x4, ν) = emissivity_absorption(obj, x4, ν)[1]
@inline absorption(obj::AbstractMedium, x4, ν) = emissivity_absorption(obj, x4, ν)[2]


# --- Synchrotron (Stage 1: angle-averaged Stokes I) ---
struct PowerLawElectrons{Tp,Tγ,TC}
	p::Tp
	γmin::Tγ
	γmax::Tγ
	Cj::TC
	Ca::TC
end

function PowerLawElectrons(; p, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)
	if isnothing(Cj) || isnothing(Ca)
		@assert isnothing(Cj) && isnothing(Ca)
		K_per_ne = _K_per_ne(p, γmin, γmax)
		(c5, c6) = _synchrotron_c5_c6(p)
		Cj = c5 * K_per_ne
		Ca = c6 * K_per_ne
	end
	return PowerLawElectrons(p, promote(γmin, γmax)..., promote(Cj, Ca)...)
end

@unstable prepare_for_computations(model::PowerLawElectrons) = @modify(FixedExponent, model.p)

# Pitch-angle average for an isotropic electron momentum distribution.
# With θ the pitch angle, isotropy implies p(θ) = (1/2)sinθ on [0, π]. Then
#   ⟨sin^q θ⟩ = (1/2)∫₀^π sin^{q+1}θ dθ
#            = √π · Γ((q+2)/2) / (2 · Γ((q+3)/2)).
@inline _avg_sin_pow(q) = sqrt(pi) * SpecialFunctions.gamma((q + 2) / 2) / (2 * SpecialFunctions.gamma((q + 3) / 2))

# Normalization for a power-law electron distribution N(γ) = K γ^{-p} on [γmin, γmax].
# The number density is
#   n_e = ∫ N(γ) dγ = K/(p-1) · (γmin^(1-p) - γmax^(1-p))  (for p>1).
# Therefore K/n_e = (p-1)/(γmin^(1-p) - γmax^(1-p)); for γmax=Inf, the γmax term → 0.
@inline _K_per_ne(p, γmin, γmax) = begin
	@assert p > 1 && γmin > 0 && γmax > γmin
	invnorm = isfinite(γmax) ? (γmin^(1 - p) - γmax^(1 - p)) : γmin^(1 - p)
	return (p - 1) / invnorm
end

# Analytic coefficient helpers (CGS); see e.g. Rybicki & Lightman.
# Returns the pair (c5(p), c6(p)) used for power-law synchrotron emissivity/absorption.
@inline _synchrotron_c5_c6(p) = begin
	e = 4.8032068e-10           # statC
	me = ustrip(u"g", u"me")     # g
	c = ustrip(u"cm/s", 1.0u"c") # cm/s

	# dimensionless conversion factor from ν/B
	A = (2 * pi * me * c) / (3 * e)
	pref0 = sqrt(3) * e^3 / (16 * pi * me)

	# pitch-angle averaged B_perp^( (p+1)/2 ) and B_perp^( (p+2)/2 )
	sinavg5 = _avg_sin_pow((p + 1) / 2)
	sinavg6 = _avg_sin_pow((p + 2) / 2)

	c5 = (pref0 / c^2) * (p - 1) * SpecialFunctions.gamma((3p - 1) / 12) * SpecialFunctions.gamma((3p + 7) / 12) * A^(-(p - 1) / 2) * sinavg5
	c6 = pref0 * (p + 2) * SpecialFunctions.gamma((3p + 2) / 12) * SpecialFunctions.gamma((3p + 10) / 12) * A^(-(p + 2) / 2) * sinavg6
	return c5, c6
end

@inline _synchrotron_coeffs(model::PowerLawElectrons, n_e, B, ν) = let
	# Stage 1 (angle-averaged) power-law synchrotron scaling, in the comoving frame.
	# This returns (j_ν, α_ν) with frequency ν measured in the plasma rest frame.
	#
	# Implemented laws:
	#   j_ν = Cj · n_e · B^((p+1)/2) · ν^(-(p-1)/2)
	#   α_ν = Ca · n_e · B^((p+2)/2) · ν^(-(p+4)/2)
	#
	# Algebraic factorization used below:
	#   common ≡ (B/ν)^(p/2)
	#   j_ν = Cj · n_e · common · √B · √ν
	#   α_ν = Ca · n_e · common · B · ν^-2
	(;p, Cj, Ca) = model
	invν = inv(ν)
	B_over_ν = B * invν
	common = B_over_ν^_half(p)
	j = Cj * n_e * common * sqrt(B*ν)
	α = Ca * n_e * common * B * invν^2
	return j, α
end

@inline emissivity_absorption(obj::AbstractSynchrotronMedium, x4, ν) =
	_synchrotron_coeffs(
		synchrotron_model(obj),
		electron_density(obj, x4), magnetic_field_strength(obj, x4), ν)

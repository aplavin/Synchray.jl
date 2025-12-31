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
@kwdef struct PowerLawElectrons{Tp,TEmin,TEmax,TCj,TCa}
	p::Tp
	Emin::TEmin = nothing
	Emax::TEmax = nothing
	Cj::TCj = 1
	Ca::TCa = 1
end

@unstable prepare_for_computations(model::PowerLawElectrons) = @modify(FixedExponent, model.p)

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

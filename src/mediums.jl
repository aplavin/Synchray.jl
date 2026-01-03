abstract type AbstractMedium end

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
@inline emissivity_absorption_invariant(obj::AbstractMedium, x4, k′) = begin
	ν′ = k′.t
	(j, α) = emissivity_absorption(obj, x4, k′)
	return (j / (ν′^2), α * ν′)
end

# if a medium defines only `emissivity_absorption`, the generic `emissivity`/`absorption` wrappers below will work
# `absorption` is used in optical depth calculation, maybe can drop this in future?..
@inline emissivity_invariant(obj::AbstractMedium, x4, k′) = emissivity(obj, x4, k′) / (k′.t^2)
@inline absorption_invariant(obj::AbstractMedium, x4, k′) = absorption(obj, x4, k′) * k′.t
@inline emissivity(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[1]
@inline absorption(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[2]



abstract type AbstractSynchrotronMedium <: AbstractMedium end

abstract type AbstractMagneticField end

"""Magnetic field that is isotropically tangled on unresolved scales (direction averaged)."""
struct FullyTangled{T<:Number} <: AbstractMagneticField
	strength::T
end
Base.:≈(a::FullyTangled, b::FullyTangled; kwargs...) = isapprox(a.strength, b.strength; kwargs...)
Base.:*(a::Number, b::FullyTangled) = FullyTangled(a * b.strength)

"""Magnetic field with a preferred direction and partial ordering parameter κ."""
struct PartiallyTangled{Tb<:SVector{3},Tκ<:Number} <: AbstractMagneticField
	b::Tb
	kappa::Tκ
end

PartiallyTangled(b; kappa) = PartiallyTangled(b, kappa)

@inline emissivity_absorption(obj::AbstractSynchrotronMedium, x4, k′) =
	_synchrotron_coeffs(
		synchrotron_model(obj),
		electron_density(obj, x4), magnetic_field(obj, x4), k′)

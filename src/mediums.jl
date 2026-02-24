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
	ν′ = frequency(k′)
	(j, α) = emissivity_absorption(obj, x4, k′)
	return (j / (ν′^2), α * ν′)
end

@inline emissivity_absorption_polarized_invariant(obj::AbstractMedium, x4, k′) = begin
	ν′ = frequency(k′)
	(j, α, B′) = emissivity_absorption_polarized(obj, x4, k′)
	return (j / (ν′^2), α * ν′, B′)
end

# if a medium defines only `emissivity_absorption`, the generic `emissivity`/`absorption` wrappers below will work
# `absorption` is used in optical depth calculation, maybe can drop this in future?..
@inline emissivity_invariant(obj::AbstractMedium, x4, k′) = emissivity(obj, x4, k′) / frequency(k′)^2
@inline absorption_invariant(obj::AbstractMedium, x4, k′) = absorption(obj, x4, k′) * frequency(k′)
@inline emissivity(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[1]
@inline absorption(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[2]



abstract type AbstractMagneticField end

"""Magnetic field that is isotropically tangled on unresolved scales (direction averaged)."""
struct FullyTangled{T<:Number} <: AbstractMagneticField
	strength::T
end
Base.:≈(a::FullyTangled, b::FullyTangled; kwargs...) = isapprox(a.strength, b.strength; kwargs...)
Base.:*(a::Number, b::FullyTangled) = FullyTangled(a * b.strength)
Base.:*(a::FullyTangled, b::Number) = FullyTangled(a.strength * b)
Base.zero(::FullyTangled) = FullyTangled(0)

"""
Magnetic field described by a mixture of:

- a fully tangled (isotropic) component, and
- a preferred ordered direction `b`.

The mixture weight is controlled by a heuristic parameter `kappa`.

This type is meant to represent “some ordering” without committing (yet) to a specific
continuous angular distribution (e.g. von Mises / Fisher with a concentration parameter).

Semantics (current minimal model):

- `b` is the preferred/ordered comoving field vector `B′`.
- `kappa ≥ 0` is *not* a concentration parameter.
- In the isotropic power-law synchrotron model (`IsotropicPowerLawElectrons`), `kappa` is mapped to a mixing fraction
	`f = kappa/(1+kappa)` (with `kappa=Inf` treated as `f=1`).
	Radiative coefficients then interpolate between:
	- `kappa = 0`: isotropically tangled (direction-averaged), equivalent to using
		`FullyTangled(|B′|)` for Stokes-I coefficients.
	- `kappa → ∞`: fully ordered, equivalent to using the raw ordered field vector `b`.
"""
struct TangledOrderedMixture{Tb<:SVector{3},Tκ<:Number} <: AbstractMagneticField
	b::Tb
	kappa::Tκ
end

TangledOrderedMixture(b; kappa) = TangledOrderedMixture(b, kappa)

Base.:≈(a::TangledOrderedMixture, b::TangledOrderedMixture; kwargs...) =
	# XXX: we propagate atol and rtol to both fields, but this don't really make sense...
	isapprox(a.b, b.b; kwargs...) && isapprox(a.kappa, b.kappa; kwargs...)

Base.:*(a::Number, f::TangledOrderedMixture) = TangledOrderedMixture(a * f.b, f.kappa)
Base.:*(f::TangledOrderedMixture, a::Number) = TangledOrderedMixture(f.b * a, f.kappa)


"""
    CombinedMedium{T<:Tuple} <: AbstractMedium

Multiple non-overlapping media rendered as a single object.
Each component gets `ray.nz` integration steps.

Objects must be ordered back-to-front (increasing z, i.e. from source toward observer).
This is the direction of photon propagation: the first object in the tuple is furthest
from the camera, and its light is attenuated by all subsequent objects.

Construct with positional arguments: `CombinedMedium(obj_back, obj_mid, obj_front)`.
"""
struct CombinedMedium{T<:Tuple} <: AbstractMedium
	objects::T
end
CombinedMedium(objs...) = CombinedMedium(objs)

@unstable prepare_for_computations(cm::CombinedMedium) = CombinedMedium(map(prepare_for_computations, cm.objects))
@unstable ustrip(cm::CombinedMedium) = CombinedMedium(map(ustrip, cm.objects))


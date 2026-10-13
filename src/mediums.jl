abstract type AbstractMedium end

# Primary API (what the raytracer should prefer calling): compute both comoving j_ν and α_ν
# (they typically share expensive intermediate computations).
#
# The integrator works in frequency-scaled transfer quantities, normalized by a per-ray
# reference frequency ν_ref (the observer frequency) so they sit at the physical-intensity
# scale and stay within Float32 range:
#   𝓘̂ ≡ I_ν·(ν_ref/ν)³   (carried; equals physical I_ν at the observer, where ν = ν_ref)
#   𝓙̂ ≡ j_ν·ν_ref³/ν²    (emission)
#   𝓐  ≡ α_ν·ν            (absorption / optical depth; independent of ν_ref)
# along the affine parameter λ:  d𝓘̂/dλ = 𝓙̂ - 𝓐·𝓘̂.
#
# Returns (𝓙̂, 𝓐) from the comoving (j_ν, α_ν) at the local frequency ν′ = frequency(k′).
# 𝓙̂ is grouped as j·ν_ref·(ν_ref/ν′)² to keep every intermediate in range — the expanded
# j·ν_ref³/ν′² under/overflows in Float32. Keep this grouping.
@inline emissivity_absorption_invariant(obj::AbstractMedium, x4, k′, ν_ref) = begin
	ν′ = frequency(k′)
	(j, α) = emissivity_absorption(obj, x4, k′)
	return (j * ν_ref * (ν_ref / ν′)^2, α * ν′)
end

@inline emissivity_absorption_polarized_invariant(obj::AbstractMedium, x4, k′, ν_ref) = begin
	ν′ = frequency(k′)
	(j, α, B′) = emissivity_absorption_polarized(obj, x4, k′)
	return (j * ν_ref * (ν_ref / ν′)^2, α * ν′, B′)
end

# if a medium defines only `emissivity_absorption`, the generic `emissivity`/`absorption` wrappers below will work
# `absorption` is used in optical depth calculation
@inline emissivity_invariant(obj::AbstractMedium, x4, k′) = emissivity(obj, x4, k′) / frequency(k′)^2
@inline absorption_invariant(obj::AbstractMedium, x4, k′) = absorption(obj, x4, k′) * frequency(k′)
@inline emissivity(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[1]
@inline absorption(obj::AbstractMedium, x4, k′) = emissivity_absorption(obj, x4, k′)[2]


"""
    is_time_independent(x) -> Bool

Trait: `true` if the object's `is_inside`, `four_velocity`, and `emissivity_absorption`
(or their leaf equivalents) do not depend on `x4.t` (the time coordinate).
Default is `true`; override to `false` for time-dependent leaves.

GR integrators may use this to skip computing the emission-time integral `t_em`,
which removes substantial work in the inner τ-loop (saves the It / Gt / I1 / I2 chain
in `SynchrayGR`).

Time-dependent leaves currently flagged: worldline-based geometries
(`Ball`, `Ellipsoid`, `WorldtubeEllipsoid`) and `Patterns.PrecessingNozzle`.
Composite types (`EmissionRegion`, `CombinedMedium`, `BFieldSpec`, ...) are
the AND of their components.
"""
@inline is_time_independent(::Any) = true



abstract type AbstractMagneticField end

"""Comoving-frame magnetic field that is isotropically tangled on unresolved scales (direction averaged)."""
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
	# atol/rtol (if given) are applied independently to both fields
	isapprox(a.b, b.b; kwargs...) && isapprox(a.kappa, b.kappa; kwargs...)

Base.:*(a::Number, f::TangledOrderedMixture) = TangledOrderedMixture(a * f.b, f.kappa)
Base.:*(f::TangledOrderedMixture, a::Number) = TangledOrderedMixture(f.b * a, f.kappa)


"""
    CombinedMedium <: AbstractMedium

Multiple non-overlapping media rendered as a single object.

Flat-space `integrate_ray` gives each intersected component `ray.nz` equal-width midpoint
integration cells. Pointwise (e.g. GR geodesic stepping) it acts as the component containing
the point; a transfer step applies the step of every component containing the point.

Objects must be ordered back-to-front (increasing z, i.e. from source toward observer).
This is the direction of photon propagation: the first object in the tuple is furthest
from the camera, and its light is attenuated by all subsequent objects.

Construct with positional arguments: `CombinedMedium(obj_back, obj_mid, obj_front)`,
or pass a vector: `CombinedMedium([obj_back, obj_mid, obj_front])`.
"""
struct CombinedMedium{T <: Union{Tuple,AbstractVector}} <: AbstractMedium
	objects::T
end
CombinedMedium(objs...) = CombinedMedium(objs)

@inline is_time_independent(cm::CombinedMedium) = all(is_time_independent, cm.objects)

@unstable prepare_for_computations(cm::CombinedMedium) = CombinedMedium(map(prepare_for_computations, cm.objects))
@unstable ustrip(cm::CombinedMedium) = CombinedMedium(map(ustrip, cm.objects))

@inline is_inside(cm::CombinedMedium, x4) = any(obj -> is_inside(obj, x4), cm.objects)

# Pointwise: the first component containing x4 (one frame cannot represent overlapping components
# moving differently). Outside all: vacuum coefficients in the first component's frame.
@inline four_velocity(cm::CombinedMedium, x4) =
	_with_inside(obj -> four_velocity(obj, x4), cm.objects, x4, () -> four_velocity(first(cm.objects), x4))

@inline emissivity_absorption(cm::CombinedMedium, x4, k′) =
	_with_inside(obj -> emissivity_absorption(obj, x4, k′), cm.objects, x4, () -> (zero(frequency(k′)), zero(frequency(k′))))

@inline emissivity_absorption_polarized(cm::CombinedMedium, x4, k′) =
	_with_inside(obj -> emissivity_absorption_polarized(obj, x4, k′), cm.objects, x4, () -> begin
		(j, α, B′) = emissivity_absorption_polarized(first(cm.objects), x4, k′)
		(zero(j), zero(α), B′)
	end)

# f(first component containing x4), or orelse() if none. Tuples unroll: no dynamic dispatch for
# heterogeneous components (GPU).
@inline _with_inside(f, ::Tuple{}, x4, orelse) = orelse()
@inline _with_inside(f, objs::Tuple, x4, orelse) =
	is_inside(first(objs), x4) ? f(first(objs)) : _with_inside(f, Base.tail(objs), x4, orelse)
@inline _with_inside(f, objs::AbstractVector, x4, orelse) = begin
	for obj in objs
		is_inside(obj, x4) && return f(obj)
	end
	orelse()
end

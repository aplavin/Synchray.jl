"""
    PeakedEmission{TS, Ta, Tν0, Tσ}

Emission model with a bell-shaped spectral profile peaked around a center frequency.
Useful for simulating colored objects that emit/absorb at specific frequencies.

The spectral shape is an algebraic bell curve, symmetric on the log-frequency axis:
    x = (ν/ν₀ - ν₀/ν) / 2    # = sinh(log(ν/ν₀))
    S(ν) = 1 / (1 + (x/σ)²)

The spectral shape multiplies the source function S_ν = j_ν/α_ν (not α_ν directly),
so that optically thick objects remain colored via Kirchhoff's law:
    α_ν = α_amp
    j_ν = α_ν · S_amp · spectral(ν)

# Fields
- `S::TS`: Source function amplitude profile, callable as `S(geom, x4) -> scalar`
- `α::Ta`: Absorption coefficient amplitude profile, callable as `α(geom, x4) -> scalar`
- `ν₀::Tν0`: Center frequency
- `σ::Tσ`: Width parameter (dimensionless)
"""
@kwdef struct PeakedEmission{TS, Ta, Tν0, Tσ}
    S::TS
    α::Ta
    ν₀::Tν0
    σ::Tσ
end

@inline function emissivity_absorption(em::PeakedEmission, geom, x4, k′)
    ν = frequency(k′)
    r = ν / em.ν₀
    x = (r - inv(r)) / 2
    spectral = inv(1 + (x / em.σ)^2)
    α_val = em.α(geom, x4)
    j_val = α_val * em.S(geom, x4) * spectral
    return (j_val, α_val)
end

@unstable prepare_for_computations(em::PeakedEmission) = modify(prepare_for_computations, em, @o _.S _.α)

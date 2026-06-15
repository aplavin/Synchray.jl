"""
    SynchrotronEmission{Tne, TB, Te}

Synchrotron emission model with electron density, magnetic field, and electron distribution.

All physical quantities (`ne`, `B`) are interpreted as **comoving-frame** (plasma rest frame)
values. The radiative transfer code handles the Lorentz transformation of the photon
frequency and the Doppler boosting of the resulting intensity automatically.

# Fields
- `ne::Tne`: Comoving electron density profile, callable as `ne(geom, x4) -> scalar`
- `B::TB`: Comoving magnetic field specification (`BFieldSpec`)
- `electrons::Te`: Electron distribution model (e.g. `IsotropicPowerLawElectrons`)
"""
@kwdef struct SynchrotronEmission{Tne, TB, Te}
    ne::Tne
    B::TB
    electrons::Te
end

@inline electron_density(em::SynchrotronEmission, geom, x4) = em.ne(geom, x4)

@inline function magnetic_field(em::SynchrotronEmission, geom, x4)
    b_mag = em.B.scale(geom, x4)
    b_dir = field_direction(em.B.dir, geom, x4)
    return em.B.wrap(b_mag * b_dir)
end

@inline emissivity_absorption(em::SynchrotronEmission, geom, x4, k′) =
    _synchrotron_coeffs(em.electrons, electron_density(em, geom, x4), magnetic_field(em, geom, x4), k′, geom, x4)

@inline emissivity_absorption_polarized(em::SynchrotronEmission, geom, x4, k′) = begin
    (jI, αI) = emissivity_absorption(em, geom, x4, k′)
    field = magnetic_field(em, geom, x4)
    (j, α) = _emissivity_absorption_polarized_field(em.electrons, jI, αI, field, k′)
    return (j, α, field)
end

# Fallback: drop (geom, x4) for models that don't need spatial context.
# BrokenPowerLaw defines a more specific 6-arg method that uses them.
@inline _synchrotron_coeffs(model, n_e, field, k′, geom, x4) = _synchrotron_coeffs(model, n_e, field, k′)

@inline is_time_independent(em::SynchrotronEmission) =
    is_time_independent(em.ne) && is_time_independent(em.B) && is_time_independent(em.electrons)

@unstable prepare_for_computations(em::SynchrotronEmission) = modify(prepare_for_computations, em, @o _.ne _.B _.electrons)
@unstable ustrip(em::SynchrotronEmission) = @p let
    em
    modify(ustrip, __, @o _.B _.electrons)
    @modify(ustrip(_; valu=UCTX.ne0), __.ne)
end

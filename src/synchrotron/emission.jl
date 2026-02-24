"""
    SynchrotronEmission{Tne, TB, Te}

Synchrotron emission model with electron density, magnetic field, and electron distribution.

# Fields
- `ne::Tne`: Electron density profile, callable as `ne(geom, x4) -> scalar`
- `B::TB`: Magnetic field specification (`BFieldSpec`)
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
    _synchrotron_coeffs(em.electrons, electron_density(em, geom, x4), magnetic_field(em, geom, x4), k′)

@inline emissivity_absorption_polarized(em::SynchrotronEmission, geom, x4, k′) = begin
    (jI, αI) = emissivity_absorption(em, geom, x4, k′)
    field = magnetic_field(em, geom, x4)
    (j, α) = _emissivity_absorption_polarized_field(em.electrons, jI, αI, field, k′)
    return (j, α, field)
end

@unstable prepare_for_computations(em::SynchrotronEmission) = modify(prepare_for_computations, em, @o _.ne _.B _.electrons)
@unstable ustrip(em::SynchrotronEmission) = @p let
    em
    modify(ustrip, __, @o _.B _.electrons)
    @modify(ustrip(_; valu=UCTX.ne0), __.ne)
end

"""
    BFieldSpec{Tscale, Tdir, Twrap}

Specification for magnetic field with direction, magnitude profile, and wrapping.

# Fields
- `scale::Tscale`: Profile that returns scalar field strength, callable as `scale(geom, x4)`
- `dir::Tdir`: Direction type from `Directions` module
- `wrap::Twrap`: Function that wraps field magnitude into `AbstractMagneticField` type

# Example
```julia
B = BFieldSpec(
    Profiles.Axial(s -> s^-1),
    Directions.Helical(π/4),
    b -> TangledOrderedMixture(b; kappa=5)
)
```
"""
struct BFieldSpec{Tscale,Tdir,Twrap}
    scale::Tscale
    dir::Tdir
    wrap::Twrap
end

prepare_for_computations(bspec::BFieldSpec) = modify(prepare_for_computations, bspec, @o _.dir _.scale _.wrap)
ustrip(bspec::BFieldSpec) = @modify(s -> ustrip(s; valu=UCTX.B0), bspec.scale)

"""
    VelocitySpec{Tdir, Tkind, Tscale}

Specification for velocity field with direction, speed type, and magnitude profile.

# Fields
- `dir::Tdir`: Direction type from `Directions` module
- `kind::Tkind`: Function specifying speed type (`beta` or `gamma`)
- `scale::Tscale`: Profile that returns scalar speed value, callable as `scale(geom, x4)`

# Constructors
```julia
VelocitySpec(dir, scale)  # defaults to gamma
VelocitySpec(dir, kind, scale)
```

# Example
```julia
# Constant Lorentz factor
velocity = VelocitySpec(Directions.Axial(), Profiles.Constant(10.0))

# Beta varying with position
velocity = VelocitySpec(
    Directions.Radial(),
    beta,
    Profiles.Natural(c -> 0.9 * (1 - c.η^2))
)
```
"""
struct VelocitySpec{Tdir, Tkind, Tscale}
    dir::Tdir
    kind::Tkind
    scale::Tscale
end

# Default constructor with gamma
VelocitySpec(dir, scale) = VelocitySpec(dir, gamma, scale)

prepare_for_computations(vspec::VelocitySpec) = modify(prepare_for_computations, vspec, @o _.dir _.kind _.scale)
ustrip(vspec::VelocitySpec) = vspec

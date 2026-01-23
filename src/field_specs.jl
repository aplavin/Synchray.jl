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
    Directions.HelicalAT(π/4),
    b -> TangledOrderedMixture(b; kappa=5)
)
```
"""
struct BFieldSpec{Tscale,Tdir,Twrap}
    scale::Tscale
    dir::Tdir
    wrap::Twrap
end

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

"""
    CombinedVelocity{T1, T2}

Sum of two velocity specifications. β-vectors are added in the lab frame,
then the 4-velocity is constructed from the combined β.

# Example
```julia
# Radial + rotation
velocity = VelocitySpec(Directions.Radial(), beta, Profiles.Transverse(β_cross)) +
           VelocitySpec(Directions.Toroidal(), beta, Profiles.RigidRotation(0.1, 1u"pc"))
```
"""
struct CombinedVelocity{T1, T2}
    v1::T1
    v2::T2
end

Base.:(+)(v1::VelocitySpec, v2::VelocitySpec) = CombinedVelocity(v1, v2)
Base.:(+)(v1::VelocitySpec, v2::CombinedVelocity) = CombinedVelocity(v1, v2)
Base.:(+)(v1::CombinedVelocity, v2::VelocitySpec) = CombinedVelocity(v1, v2)
Base.:(+)(v1::CombinedVelocity, v2::CombinedVelocity) = CombinedVelocity(v1, v2)


@unstable begin
prepare_for_computations(bspec::BFieldSpec) = modify(prepare_for_computations, bspec, @o _.dir _.scale _.wrap)
prepare_for_computations(vspec::VelocitySpec) = modify(prepare_for_computations, vspec, @o _.dir _.kind _.scale)
prepare_for_computations(vel::CombinedVelocity) = CombinedVelocity(prepare_for_computations(vel.v1), prepare_for_computations(vel.v2))
ustrip(bspec::BFieldSpec) = @modify(s -> ustrip(s; valu=UCTX.B0), bspec.scale)
ustrip(vspec::VelocitySpec) = @modify(s -> ustrip(s; valu=1), vspec.scale)
ustrip(vel::CombinedVelocity) = CombinedVelocity(ustrip(vel.v1), ustrip(vel.v2))
end

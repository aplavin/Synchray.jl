"""
Module for profile wrapper types.

Profiles define how physical quantities (density, field strength, velocity)
vary spatially, with different patterns of dependence on natural coordinates.
"""
module Profiles

# No type exports - always use Profiles.Axial(), etc.

"""
    Axial{F}

Profile depending only on the axial coordinate `s`.

# Example
```julia
ne = Profiles.Axial(PowerLaw(-2; val0=1e4, s0=1pc))
```
"""
struct Axial{F}
    f::F
end

"""
    Transverse{F}

Profile depending only on the normalized transverse coordinate `η`.

# Example
```julia
ne = Profiles.Transverse(η -> exp(-η^2))
```
"""
struct Transverse{F}
    f::F
end

"""
    AxialTransverse{Fz, Fη}

Profile as a separable product of axial and transverse functions: `f(z) * g(η)`.

# Example
```julia
ne = Profiles.AxialTransverse(
    z -> (z/z0)^-2,           # axial decline
    η -> exp(-η^2 / 0.5^2)    # Gaussian transverse profile
)
```
"""
struct AxialTransverse{Fz, Fη}
    f_z::Fz      # function of z
    f_η::Fη      # function of η
end

"""
    Natural{F}

Profile depending on the full natural coordinate tuple.

# Example
```julia
ne = Profiles.Natural(c -> c.z^-2 * exp(-c.η^2))
```
"""
struct Natural{F}
    f::F
end

"""
    Raw{F}

Profile with raw access to geometry and position.

# Example
```julia
ne = Profiles.Raw((geom, x4) -> custom_density(geom, x4))
```
"""
struct Raw{F}
    f::F
end

"""
    Constant{T}

Constant value profile (independent of position).

# Example
```julia
ne = Profiles.Constant(1e4)
```
"""
struct Constant{T}
    val::T
end

"""
    Modified{Tbase, Tmod}

Profile modified by a multiplicative pattern/modifier.

# Fields
- `base`: Original profile
- `modifier`: Function `(geom, x4, base_val) -> modified_val`

# Example
```julia
ne = Profiles.Modified(
    Profiles.Axial(base_density),
    (geom, x4, val) -> val * pattern_factor(geom, x4)
)
```
"""
struct Modified{Tbase, Tmod}
    base::Tbase
    modifier::Tmod
end

end # module Profiles

@inline function (p::Profiles.Axial)(geom, x4)
    z = natural_coords(geom, x4, Val(:z))
    return p.f(z)
end
prepare_for_computations(p::Profiles.Axial) = modify(prepare_for_computations, p, @o _.f)
ustrip(p::Profiles.Axial; valu) = @modify(f -> ustrip(f; valu, argu=UCTX.L0), p.f)

@inline function (p::Profiles.Transverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f(coords.η)
end
prepare_for_computations(p::Profiles.Transverse) = modify(prepare_for_computations, p, @o _.f)
ustrip(p::Profiles.Transverse; valu) = @modify(f -> ustrip(f; valu, argu=1), p.f)

@inline function (p::Profiles.AxialTransverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f_z(coords.z) * p.f_η(coords.η)
end
prepare_for_computations(p::Profiles.AxialTransverse) = modify(prepare_for_computations, p, @o _.f_z _.f_η)
ustrip(p::Profiles.AxialTransverse; valu) = @p let
    p
    @modify(ustrip(_; valu, argu=UCTX.L0), __.f_z)
    @modify(ustrip(_; valu, argu=1), __.f_η)
end

@inline function (p::Profiles.Natural)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f(coords)
end

@inline (p::Profiles.Raw)(geom, x4) = p.f(geom, x4)

@inline (p::Profiles.Constant)(geom, x4) = p.val

@inline function (p::Profiles.Modified)(geom, x4)
    base_val = p.base(geom, x4)
    return p.modifier(geom, x4, base_val)
end
prepare_for_computations(p::Profiles.Modified) = modify(prepare_for_computations, p, @o _.base _.modifier)
ustrip(p::Profiles.Modified; valu) = @p let
    p
    @modify(ustrip(_; valu), __.base)
end

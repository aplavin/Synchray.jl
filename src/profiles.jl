"""
Module for profile wrapper types.

Profiles define how physical quantities (density, field strength, velocity)
vary spatially, with different patterns of dependence on natural coordinates.
"""
module Profiles

# Will be available from Geometries module after both are loaded
using ..Geometries: natural_coords

# No type exports - always use Profiles.Axial(), etc.

"""
    Axial{F}

Profile depending only on the axial coordinate `s`.

# Example
```julia
ne = Profiles.Axial(PowerLawS(-2; val0=1e4, s0=1pc))
```
"""
struct Axial{F}
    f::F
end

function (p::Axial)(geom, x4)
    s = natural_coords(geom, x4, Val(:s))
    return p.f(s)
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

function (p::Transverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f(coords.η)
end

"""
    AxialTransverse{Fs, Fη}

Profile as a separable product of axial and transverse functions: `f(s) * g(η)`.

# Example
```julia
ne = Profiles.AxialTransverse(
    s -> (s/s0)^-2,           # axial decline
    η -> exp(-η^2 / 0.5^2)    # Gaussian transverse profile
)
```
"""
struct AxialTransverse{Fs, Fη}
    f_s::Fs      # function of s
    f_η::Fη      # function of η
end

function (p::AxialTransverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f_s(coords.s) * p.f_η(coords.η)
end

"""
    Natural{F}

Profile depending on the full natural coordinate tuple.

# Example
```julia
ne = Profiles.Natural(c -> c.s^-2 * exp(-c.η^2))
```
"""
struct Natural{F}
    f::F
end

function (p::Natural)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f(coords)
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

(p::Raw)(geom, x4) = p.f(geom, x4)

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

(p::Constant)(geom, x4) = p.val

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

function (p::Modified)(geom, x4)
    base_val = p.base(geom, x4)
    return p.modifier(geom, x4, base_val)
end

end # module Profiles

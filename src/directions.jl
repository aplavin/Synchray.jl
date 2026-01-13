"""
Module for field direction types.

Direction types specify how field vectors (B-field, velocity) are oriented 
in space relative to the geometry.
"""
module Directions

using StaticArrays
using LinearAlgebra
using ..Geometries

# Export only interface functions, not type names
export field_direction

"""
    AbstractDirection

Abstract type for field direction specifications.
"""
abstract type AbstractDirection end

"""
    Scalar <: AbstractDirection

No directional component - returns scalar 1.
Used for isotropically tangled fields where only the magnitude matters.
"""
struct Scalar <: AbstractDirection end

"""
    Axial <: AbstractDirection

Field direction along the geometry's primary axis.
"""
struct Axial <: AbstractDirection end

"""
    Radial <: AbstractDirection

Field direction radially outward from the origin.
"""
struct Radial <: AbstractDirection end

"""
    Toroidal <: AbstractDirection

Field direction perpendicular to both the axis and radial direction.
Forms azimuthal loops around the axis.
"""
struct Toroidal <: AbstractDirection end

"""
    Helical{T} <: AbstractDirection

Field direction as a helical mix of axial and toroidal components.

# Fields
- `ψ::T`: Pitch angle (angle between field and axis)
"""
struct Helical{T} <: AbstractDirection
    ψ::T  # pitch angle
end

"""
    field_direction(dir::AbstractDirection, geometry, x4) -> SVector{3} or scalar

Compute the field direction unit vector (or scalar for Scalar direction) 
at spacetime position `x4` for the given geometry.

For `Scalar`, returns `1` (no direction).
For vector directions, returns a unit 3-vector in the lab frame.
"""
function field_direction end

# Scalar case: return 1
@inline field_direction(::Scalar, geom, x4) = 1

# Axial: along geometry axis
@inline field_direction(::Axial, geom, x4) = Geometries.geometry_axis(geom)

# Radial: outward from origin (normalized position vector)
@inline field_direction(::Radial, geom, x4) = begin
    r = SVector(x4.x, x4.y, x4.z)
    r_norm = norm(r)
    return r / r_norm
end

# Toroidal: azimuthal around axis
@inline field_direction(::Toroidal, geom, x4) = begin
    axis = Geometries.geometry_axis(geom)
    r = SVector(x4.x, x4.y, x4.z)
    s = dot(axis, r)
    r_perp = r - s * axis
    ρ = norm(r_perp)
    
    if iszero(ρ)
        return zero(axis)
    end
    
    e_r = r_perp / ρ
    return cross(axis, e_r)
end

# Helical: mix of axial and toroidal
@inline field_direction(h::Helical, geom, x4) = begin
    axis = Geometries.geometry_axis(geom)
    r = SVector(x4.x, x4.y, x4.z)
    s = dot(axis, r)
    r_perp = r - s * axis
    ρ = norm(r_perp)
    
    e_phi = if iszero(ρ)
        zero(axis)
    else
        cross(axis, r_perp / ρ)
    end
    
    sψ, cψ = sincos(h.ψ)
    v = cψ * axis + sψ * e_phi
    return iszero(v) ? zero(v) : v / norm(v)
end

end # module Directions

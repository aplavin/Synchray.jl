"""
Module for field direction types.

Direction types specify how field vectors (B-field, velocity) are oriented 
in space relative to the geometry.
"""
module Directions

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

end # module Directions

"""
    field_direction(dir::AbstractDirection, geometry, x4) -> SVector{3} or scalar

Compute the field direction unit vector (or scalar for Scalar direction) 
at spacetime position `x4` for the given geometry.

For `Scalar`, returns `1` (no direction).
For vector directions, returns a unit 3-vector in the lab frame.
"""
function field_direction end

# Scalar case: return 1
@inline field_direction(::Directions.Scalar, geom, x4) = 1

# Axial: along geometry axis
@inline field_direction(::Directions.Axial, geom, x4) = geometry_axis(geom)

# Radial: outward from origin (normalized position vector)
@inline field_direction(::Directions.Radial, geom, x4) = begin
    r = @swiz x4.xyz
    r_norm = norm(r)
    return iszero(r_norm) ? zero(r) : r / r_norm
end

# Toroidal: azimuthal around axis
@inline field_direction(::Directions.Toroidal, geom, x4) = begin
    axis = geometry_axis(geom)
    (; r_perp, ρ) = _cylindrical_coords(rotation_local_to_lab(geom), x4)

    return iszero(ρ) ? zero(axis) : cross(axis, r_perp / ρ)
end

prepare_for_computations(h::Directions.Helical) = @modify(Geometries.AngleTrigCached_fromangle, h.ψ)

# Helical: mix of axial and toroidal
@inline field_direction(h::Directions.Helical, geom, x4) = begin
    axis = geometry_axis(geom)
    (; r_perp, ρ) = _cylindrical_coords(rotation_local_to_lab(geom), x4)

    e_phi = iszero(ρ) ? zero(axis) : cross(axis, r_perp / ρ)

    sψ, cψ = sincos(h.ψ)
    v = cψ * axis + sψ * e_phi
    return iszero(v) ? zero(v) : v / norm(v)
end

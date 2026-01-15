"""
Module for profile wrapper types.

Profiles define how physical quantities (density, field strength, velocity)
vary spatially, with different patterns of dependence on natural coordinates.
"""
module Profiles

# No type exports - always use Profiles.Axial(), etc.

struct PowerLaw{Texp,Tval,Ts0}
	exp::Texp
	val0::Tval
	s0::Ts0
end
PowerLaw(exp; val0, s0=one(val0)) = PowerLaw(exp, val0, s0)

using DispatchDoctor: @unstable
@unstable struct LinearInterp{Tpts}
	points::Tpts  # Tuple of (x, y) pairs, sorted by x
	
    function LinearInterp(points)
        @assert !isempty(points)
		# Sort points by x coordinate
		sorted_pts = sort(collect(points); by=first) |> Tuple
		new{typeof(sorted_pts)}(sorted_pts)
	end
end

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
    Radial{F}

Profile depending only on the distance from the origin `r = |x⃗|`.

# Example
```julia
ne = Profiles.Radial(r -> exp(-r^2))
```
"""
struct Radial{F}
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


@inline (pl::Profiles.PowerLaw)(s) = s > 0 ? pl.val0 * (s / pl.s0)^pl.exp : zero(float(pl.val0))

@inline function (li::Profiles.LinearInterp)(s)
	pts = li.points
	length(pts) == 1 && return last(first(pts))
	
	x1, y1 = first(pts)
	s <= x1 && return y1  # Flat before first point
	
	xn, yn = last(pts)
	s >= xn && return yn  # Flat after last point
	
	# Linear interpolation between points
	for i in 1:length(pts)-1
		x_lo, y_lo = pts[i]
		x_hi, y_hi = pts[i+1]
		if x_lo <= s <= x_hi
			t = (s - x_lo) / (x_hi - x_lo)
			return y_lo + t * (y_hi - y_lo)
		end
	end
	return yn  # Fallback
end


@inline function (p::Profiles.Axial)(geom, x4)
    z = natural_coords(geom, x4, Val(:z))
    return p.f(z)
end

@inline function (p::Profiles.Transverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f(coords.η)
end

@inline function (p::Profiles.Radial)(geom, x4)
    r = norm(@swiz x4.xyz)
    return p.f(r)
end

@inline function (p::Profiles.AxialTransverse)(geom, x4)
    coords = natural_coords(geom, x4)
    return p.f_z(coords.z) * p.f_η(coords.η)
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


@unstable begin

prepare_for_computations(pl::Profiles.PowerLaw) = @modify(FixedExponent, pl.exp)
prepare_for_computations(p::Profiles.Modified) = modify(prepare_for_computations, p, @o _.base _.modifier)
prepare_for_computations(p::Profiles.Axial) = modify(prepare_for_computations, p, @o _.f)
prepare_for_computations(p::Profiles.Transverse) = modify(prepare_for_computations, p, @o _.f)
prepare_for_computations(p::Profiles.Radial) = modify(prepare_for_computations, p, @o _.f)
prepare_for_computations(p::Profiles.AxialTransverse) = modify(prepare_for_computations, p, @o _.f_z _.f_η)

ustrip(pl::Profiles.PowerLaw; argu, valu) = Profiles.PowerLaw(pl.exp; val0=_u_to_code(pl.val0, valu), s0=_u_to_code(pl.s0, argu))
ustrip(li::Profiles.LinearInterp; argu, valu) = Profiles.LinearInterp(map(((x, y),) -> (_u_to_code(x, argu), _u_to_code(y, valu)), li.points))
ustrip(p::Profiles.Axial; valu) = @modify(f -> ustrip(f; valu, argu=UCTX.L0), p.f)
ustrip(p::Profiles.Transverse; valu) = @modify(f -> ustrip(f; valu, argu=1), p.f)
ustrip(p::Profiles.Radial; valu) = @modify(f -> ustrip(f; valu, argu=UCTX.L0), p.f)
ustrip(p::Profiles.AxialTransverse; valu) = @p let
    p
    @modify(ustrip(_; valu, argu=UCTX.L0), __.f_z)
    @modify(ustrip(_; valu, argu=1), __.f_η)
end
ustrip(p::Profiles.Modified; valu) = @p let
    p
    @modify(ustrip(_; valu), __.base)
    @modify(ustrip, __.modifier)
end

end

"""
Module for pattern modifiers (knots, etc.) that work with EmissionRegion.

Patterns are profile transformations that modify base values multiplicatively.
They compose naturally with the Profiles.Modified wrapper.
"""
module Patterns

# No type exports - always use Patterns.EllipsoidalKnot(), etc.

"""
    FixedSizing{Ta}

Explicit fixed rest-frame semi-axes `(a_perp, a_parallel)`.

# Fields
- `a_perp::Ta`: Perpendicular semi-axis
- `a_parallel::Ta`: Parallel semi-axis

# Notes
The fields are stored as `(a_perp, a_parallel)`.
"""
struct FixedSizing{Ta}
	a_perp::Ta
	a_parallel::Ta
end

"""
    CrossSectionSizing{Tf}

Cross-section-tied sizing for conical geometries.

For conical jets the local radius scales as `R(z) = z * tan(φj)`.
This sizing chooses:
- `a_perp = f_perp * R(z_c)` - fraction of local jet radius
- `a_parallel = q * a_perp` - aspect ratio

# Fields
- `f_perp::Tf`: Fraction of local jet cross-section
- `q::Tf`: Aspect ratio `a_parallel / a_perp`

# Notes
The knot occupies a fixed fraction of the local jet cross-section.
"""
struct CrossSectionSizing{Tf}
	f_perp::Tf
	q::Tf
end

"""
    GaussianBump{Tf}

A smooth bump profile with `f(χ=0) == f_peak` and `f(χ→∞) → 1`.

# Fields
- `f_peak::Tf`: Peak multiplicative factor at χ=0
- `χ_threshold::Tf`: Threshold beyond which factor returns to 1 (default: 16)

# Physical meaning
- If used for `ne`, represents a localized overdensity/underdensity
- If used for `B`, represents a localized magnetic-field enhancement
- The profile is `f(χ) = 1 + (f_peak - 1) * exp(-χ/2)` for χ < χ_threshold
"""
@kwdef struct GaussianBump{Tf}
	f_peak::Tf
	χ_threshold::Tf = 16.0
end
GaussianBump(f_peak, χ_threshold) = GaussianBump(promote(f_peak, χ_threshold)...)
GaussianBump(f_peak) = GaussianBump(; f_peak=float(f_peak))

# In the knot frame, χ is an ellipsoidal "radius-squared" and this is a Gaussian in χ.
@inline (g::GaussianBump)(χ) =
	χ < g.χ_threshold ?
		1 + (g.f_peak - 1) * exp(-χ / 2) :
		one(float(g.f_peak))

"""
    EllipsoidalKnot{TX, TU, TS, TP}

An SR-consistent "world-tube" pattern with an ellipsoidal profile in the knot rest frame.

The knot center follows an inertial worldline with constant 4-velocity `u`:
`x_c(τ) = x_c0 + u * τ`.

At an event `x4`, we compute the knot proper time at simultaneity in the knot frame
and the spatial offset in that frame, then evaluate an axisymmetric ellipsoid.

# Fields
- `x_c0::TX`: Center position/event (FourPosition)
- `u::TU`: 4-velocity (FourVelocity)
- `sizing::TS`: Knot size parameters (FixedSizing or CrossSectionSizing)
- `profile::TP`: χ → factor function (e.g., GaussianBump)

# Usage
Use as a modifier in `Profiles.Modified`:
```julia
ne = Profiles.Modified(
    Profiles.Axial(base_density),
    Patterns.EllipsoidalKnot(
        x_c0 = FourPosition(t0, x0, y0, z0),
        u = construct(FourVelocity, gamma => 10, direction => SA[0,0,1]),
        sizing = Patterns.CrossSectionSizing(f_perp=0.3, q=2.0),
        profile = Patterns.GaussianBump(f_peak=3.0)
    )
)
```
"""
@kwdef struct EllipsoidalKnot{TX, TU, TS, TP}
	x_c0::TX
	u::TU
	sizing::TS
	profile::TP
end

end # module Patterns


# ============================================================================
# Implementation: Knot evaluation
# ============================================================================

"""Evaluate knot as a modifier: (geom, x4, base_val) -> modified_val"""
@inline function (knot::Patterns.EllipsoidalKnot)(geom, x4, base_val)
	χ = _knot_chi(knot, geom, x4)
	factor = knot.profile(χ)
	return base_val * factor
end

"""
Compute the ellipsoidal "radius-squared" χ(x) in the knot rest frame.

Steps:
1. Find the knot proper time τ such that x_c = x_c0 + u τ is simultaneous with x4
   in the knot rest frame (so u⋅(x4 - x_c) = 0).
2. Decompose Δ_rest into components parallel/transverse to the knot principal axis
   (taken to be the direction of motion) and evaluate the axisymmetric ellipsoidal
   quadratic form.
"""
@inline _knot_chi(knot::Patterns.EllipsoidalKnot, geom, x4::FourPosition) = begin
	Δ0 = x4 - knot.x_c0  # displacement from reference center event x_c0 to x4
	τ = -minkowski_dot(knot.u, Δ0)  # simultaneity proper time τ*(x4) = -u⋅(x4 - x_c0)
	x_c = knot.x_c0 + knot.u * τ  # center event x_c(τ) simultaneous with x4 in knot frame
	Δ_rest = x4 - x_c  # displacement from the simultaneous center event to x4 (u⋅Δc ≈ 0)
	# should be orthogonal to u

	e_par = _knot_velocity_axis(knot.u)  # unit principal axis in the knot rest space

	Δ_par = minkowski_dot(e_par, Δ_rest)  # component along the principal axis
	Δ_rest2 = minkowski_dot(Δ_rest, Δ_rest)  # squared rest-space distance (Δ_rest is spacelike)
	# Numerical guard: small negative values can appear from roundoff.
	Δ_perp2 = max(0, Δ_rest2 - Δ_par^2)  # squared distance transverse to e_par

	(a_par, a_perp) = _knot_sizes(knot.sizing, τ, x_c, geom)  # rest-frame semi-axes (parallel, perpendicular)
	return Δ_par^2 / a_par^2 + Δ_perp2 / a_perp^2
end

"""
Unit spacelike 4-vector aligned with the knot's lab-frame velocity direction.

Math: look for a unit spacelike vector e_par = (a, b β̂) whose spatial direction
is along β = beta(u) (β̂ = β/|β|), and which lives in the knot rest space:
u⋅e_par = 0 and e_par⋅e_par = 1 (signature (-,+,+,+)).
"""
@inline _knot_velocity_axis(u::FourVelocity) = begin
	β = beta(u)
	βmag = √dot(β, β)
	βhat = β / βmag
	γ = u.t
	return FourPosition(γ * βmag, γ * βhat)
end

"""Return sizes in the order (a_par, a_perp) to match the variables used in `_knot_chi`."""
@inline _knot_sizes(sizing::Patterns.FixedSizing, τ, x_c, geom) = (sizing.a_parallel, sizing.a_perp)

"""
Compute sizes for CrossSectionSizing with Conical geometry.

As the knot moves outward, it expands to keep a fixed fraction of the local jet radius.
"""
@inline _knot_sizes(sizing::Patterns.CrossSectionSizing, τ, x_c, geom::Geometries.Conical) = begin
	# Get axial coordinate of the current center position
	z_c = natural_coords(geom, x_c, Val(:z))
	
	# Compute rest-frame semi-axes from the local jet radius
	# Intuitively: as the knot moves outward, it expands to keep a fixed fraction of the local jet radius.
	a_perp = sizing.f_perp * (z_c * tan(geom.φj))
	a_parallel = sizing.q * a_perp
	
	return (a_parallel, a_perp)
end

# ============================================================================
# Validation and preparation
# ============================================================================

"""
Causality guard for CrossSectionSizing with Conical geometry.

With a_perp(τ) = f_perp * R(z_c(τ)) = f_perp * z_c(τ) * tan(φj), we have:
  v_perp = d a_perp / dτ = f_perp * tan(φj) * d z_c / dτ
  v_parallel = d a_parallel / dτ = q * v_perp

For inertial motion, the spatial part of a 4-velocity is u⃗ = d r⃗ / dτ, so
d z_c / dτ = geom.axis ⋅ u⃗.

Physical meaning: if the knot size is tied to the jet radius, the knot must not
expand faster than light in its own rest frame.
"""
function _validate_knot_causality(sizing::Patterns.CrossSectionSizing, u::FourVelocity, geom::Geometries.Conical)
	ds_dtau = dot(geom.axis, @swiz u.xyz)
	v_perp = abs(sizing.f_perp * tan(geom.φj) * ds_dtau)
	v_parallel = sizing.q * v_perp
	(v_perp < 1 && v_parallel < 1) ||
		error("Unphysical knot sizing: derived rest-frame expansion speeds (v_perp=$v_perp, v_par=$v_parallel) must be < 1")
	return nothing
end

"""
    validate_pattern(knot::EllipsoidalKnot, geom)

Validate that a knot pattern is physically allowed for the given geometry.

For `CrossSectionSizing` with `Conical` geometry, checks that the knot expansion
speeds in its rest frame are subluminal.

# Example
```julia
knot = Patterns.EllipsoidalKnot(...)
region = EmissionRegion(geometry=Geometries.Conical(...), ...)
validate_pattern(knot, region.geometry)  # Throws error if invalid
```
"""
function validate_pattern(knot::Patterns.EllipsoidalKnot, geom::Geometries.Conical)
	knot.sizing isa Patterns.CrossSectionSizing && _validate_knot_causality(knot.sizing, knot.u, geom)
	return nothing
end

# Default: no validation needed for other geometry types
validate_pattern(knot::Patterns.EllipsoidalKnot, geom) = nothing

"""Prepare knot for computations (propagates to fields)"""
prepare_for_computations(knot::Patterns.EllipsoidalKnot) = knot

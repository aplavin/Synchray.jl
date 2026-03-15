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
    TophatBump{Tf}

A tophat/step function profile: constant enhancement inside, base value outside.

# Fields
- `f_peak::Tf`: Multiplicative factor for χ < 1

# Physical meaning
- Used for sharp boundaries (e.g., beam edges)
- The profile is `f(χ) = f_peak` for χ < 1, otherwise `f(χ) = 1`
"""
struct TophatBump{Tf}
	f_peak::Tf
end

@inline (t::TophatBump)(χ) = χ < 1 ? t.f_peak : one(t.f_peak)

"""
    InvertedProfile{T}

Wraps a profile to compute exact pointwise reciprocals.

For any profile `p`, `InvertedProfile(p)(χ) == 1 / p(χ)` exactly.
Use `inv(profile)` to create inverted profiles.

# Example
```julia
g = GaussianBump(f_peak=3.0)
g_inv = inv(g)  # Returns InvertedProfile(g)
@assert g(χ) * g_inv(χ) ≈ 1.0  # Exact inverse at all χ
```
"""
struct InvertedProfile{T}
	inner::T
end

@inline (p::InvertedProfile)(χ) = inv(p.inner(χ))

# Inverse operations: exact pointwise reciprocals
Base.inv(g::GaussianBump) = InvertedProfile(g)
Base.inv(t::TophatBump) = TophatBump(inv(t.f_peak))  # Exact via parameter inversion
Base.inv(p::InvertedProfile) = p.inner  # Double inversion recovers original

"""
    ComplementedProfile{T}

Wraps a profile to compute its spatial complement (swaps inside/outside regions).

For bump profiles `p(χ)` with `p(0) = f_peak` and `p(∞) = 1`, the complement has
`c(0) = 1` and `c(∞) = 1/f_peak`. This is achieved via the formula: `c(χ) = p(χ) / f_peak`.

Physical meaning:
- Original: "Enhance pattern region by factor f_peak, normal elsewhere"
- Complement: "Normal in pattern region, suppress elsewhere by factor 1/f_peak"

Use `Patterns.complement(profile)` to create complemented profiles.

# Example
```julia
# TophatBump: enhance inside, normal outside
bump = TophatBump(f_peak=3.0)
@assert bump(0.5) ≈ 3.0  # Inside: enhanced
@assert bump(1.5) ≈ 1.0  # Outside: normal

# Complement: normal inside, suppress outside
comp = Patterns.complement(bump)
@assert comp(0.5) ≈ 1.0    # Inside: normal
@assert comp(1.5) ≈ 1/3    # Outside: suppressed

# Double complement recovers original
@assert Patterns.complement(Patterns.complement(bump))(χ) ≈ bump(χ)  # for all χ
```
"""
struct ComplementedProfile{T}
	inner::T
end

@inline (p::ComplementedProfile)(χ) = p.inner(χ) / p.inner.f_peak

# Complement operations
"""
    complement(profile)

Create the spatial complement of a bump profile.

For bumps with `p(0) = f_peak` and `p(∞) = 1`, returns a complemented profile with
`c(0) = 1` and `c(∞) = 1/f_peak`, effectively swapping which spatial region is modified.

# Example
```julia
bump = TophatBump(10.0)
comp = Patterns.complement(bump)
comp(0.5)  # Returns 1.0 (normal inside)
comp(1.5)  # Returns 0.1 (suppressed outside)
```
"""
complement(p) = ComplementedProfile(p)
complement(p::ComplementedProfile) = p.inner  # Double complement = identity

"""
    PrecessingNozzle{Tθ, Tφ, TT, T, TP}

A precessing nozzle pattern that creates an enhanced-emission channel rotating within a conical geometry.

The nozzle axis precesses around the geometry's primary axis, and only plasma emitted when the
nozzle was pointing near the current radial direction receives enhancement.

# Time convention
All times are in the **lab frame**. The pattern accounts for ballistic plasma flow: plasma at position `r`
observed at lab time `t` was **ejected** at `t_ej = t - r/β_flow` (flow travel time). The nozzle orientation 
at ejection time determines whether the plasma receives enhancement (frozen-flow picture).

# Fields
- `θ_precession::Tθ`: Precession angle from geometry axis (radians)
- `θ_nozzle::Tφ`: Nozzle cone half-opening angle (radians)
- `period::TT`: Precession period in lab frame code time units
- `φ0::T`: Initial phase (default: 0)
- `β_flow::T`: Flow speed for emission-time calculation (c=1 units)
- `profile::TP`: Profile function mapping normalized angle to enhancement factor

# Usage
Use as a modifier in `Profiles.Modified`:
```julia
ne = Profiles.Modified(
    Profiles.Axial(base_density),
    Patterns.PrecessingNozzle(
        θ_precession = 0.02,
        θ_nozzle = 0.005,
        period = 50.0,
        β_flow = 0.99,
        profile = Patterns.TophatBump(10.0)
    )
)
```
"""
struct PrecessingNozzle{Tθ, Tφ, TT, T, TP}
	θ_precession::Tθ
	θ_nozzle::Tφ
	period::TT
	φ0::T
	β_flow::T
	profile::TP
end
PrecessingNozzle(θ_precession, θ_nozzle, period, φ0, β_flow, profile) = PrecessingNozzle(θ_precession, θ_nozzle, period, promote(φ0, β_flow)..., profile)
PrecessingNozzle(; θ_precession, θ_nozzle, period, φ0=0, β_flow, profile) = PrecessingNozzle(θ_precession, θ_nozzle, period, φ0, β_flow, profile)

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
        u = construct(FourVelocity, gamma => 10, direction3 => SA[0,0,1]),
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

"""Evaluate precessing nozzle as a modifier: (geom, x4, base_val) -> modified_val"""
@inline function (nozzle::Patterns.PrecessingNozzle)(geom, x4, base_val)
	χ = _nozzle_chi(nozzle, geom, x4)
	factor = nozzle.profile(χ)
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
# Implementation: Rotating beam evaluation
# ============================================================================

"""
Compute the normalized angular distance from the nozzle axis for a precessing nozzle pattern.

Returns χ = sin(angle) / sin(θ_nozzle), so χ < 1 means inside the nozzle cone.
For small angles, this is approximately χ ≈ angle / θ_nozzle.
"""
@inline _nozzle_chi(nozzle::Patterns.PrecessingNozzle, geom, x4::FourPosition) = begin
	r_vec = @swiz x4.xyz
	r = norm(r_vec)
	
	# Ejection time (lab frame): when was this plasma ejected from the nozzle?
	# Plasma travels ballistically at speed β_flow, so travel time is r/β_flow
	t_ej = x4.t - r / nozzle.β_flow
	
	# Nozzle axis direction at ejection time
	nozzle_dir = _nozzle_axis_at_time(nozzle, geom, t_ej)
	
	# Angular distance from nozzle axis
	cos_angle = dot(r_vec, nozzle_dir) / r

	# For small angles, use sin(angle) / sin(θ_nozzle) ≈ angle / θ_nozzle
	# Compute sin directly from cos to avoid acos() call: sin²θ + cos²θ = 1
	sin_angle = √(max(0, 1 - cos_angle^2))
	
	# Normalize: χ = sin(angle) / sin(θ_nozzle), so χ = 0 on axis, χ = 1 at nozzle edge
	return sin_angle / sin(nozzle.θ_nozzle)
end

"""
Compute the nozzle axis direction at a given lab frame time.

The nozzle axis precesses around the geometry's primary axis with the given period.
Time `t` is in the lab frame.
"""
@inline _nozzle_axis_at_time(nozzle::Patterns.PrecessingNozzle, geom, t) = begin
	# Precession phase at this lab frame time
	phase = t * 2*π / nozzle.period + nozzle.φ0
	
	# Nozzle direction: rotate around axis by angle θ_precession
	sθ, cθ = sincos(nozzle.θ_precession)
	sφ, cφ = sincos(phase)

	# Express in local frame, then rotate to lab frame
	return rotate_local_to_lab(geom, SVector(sθ * cφ, sθ * sφ, cθ))
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
	axis = geometry_axis(geom)
	ds_dtau = dot(axis, @swiz u.xyz)
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

# Design note: For non-inertial worldlines, a cleaner design would be:
#   worldline(tau) -> FourPosition       # worldline parameterized by proper time
#   retarded_time(worldline; t_obs) -> tau   # find retarded proper time
# Current design bundles both for convenience with inertial worldlines.
"""
    retarded_event(wl::InertialWorldline, cam::CameraOrtho)
    retarded_event(knot::EllipsoidalKnot, cam::CameraOrtho)

Compute the retarded event on an inertial worldline `x(τ) = wl.x0 + wl.u·τ`
as seen by the camera `cam`.

Finds the event whose light reaches the camera's screen plane at observer time `cam.t`.
The camera defines a null hyperplane `t - n̂·r = cam.t - n̂·cam.origin`, and this function
solves for the unique intersection of that hyperplane with the worldline.

# Returns
`(; x::FourPosition, tau)` where:
- `x`: The retarded spacetime event on the worldline.
- `tau`: The retarded proper time on the worldline.
"""
@inline function retarded_event(wl::Geometries.InertialWorldline, cam::CameraOrtho)
    # Inertial worldline: x(τ) = x0 + u·τ
    # Null hyperplane: t - n̂·r = t_obs - n̂·origin
    #
    # Solve for τ:
    #   (x0.t + u.t·τ) - n̂·(x0_xyz + u_xyz·τ) = t_obs - n̂·origin
    #   τ·(u.t - n̂·u_xyz) = t_obs - n̂·origin - (x0.t - n̂·x0_xyz)
    (; x0, u) = wl
    n̂ = cam.n
    x0_xyz = @swiz x0.xyz
    u_xyz = @swiz u.xyz
    tau = (cam.t - dot(n̂, cam.origin) - (x0.t - dot(n̂, x0_xyz))) / (u.t - dot(n̂, u_xyz))
    x = x0 + u * tau
    return (; x, tau)
end

@inline retarded_event(knot::Patterns.EllipsoidalKnot, cam::CameraOrtho) =
    retarded_event(Geometries.InertialWorldline(knot.x_c0, knot.u), cam)

ustrip(knot::Patterns.EllipsoidalKnot) = @modify(x -> _u_to_code(x, UCTX.L0), knot.x_c0)

"""
Prepare precessing nozzle for computations by caching trig values.
"""
prepare_for_computations(nozzle::Patterns.PrecessingNozzle) = modify(Geometries.AngleTrigCached_fromangle, nozzle, @o _.θ_precession _.θ_nozzle)

ustrip(nozzle::Patterns.PrecessingNozzle) = modify(x -> _u_to_code(x, UCTX.T0), nozzle, @o _.period)

"""
		AbstractJetPattern

A *pattern* is a purely local, time-dependent modifier of jet microphysics.

Patterns are designed to be “slow-light correct” automatically because they are evaluated
at events `x4` that already lie on the ray worldline, with `x4.t` being **lab time**.

Contract:

- Patterns should usually modify only `electron_density` and/or `magnetic_field_strength`
	(both are comoving/proper quantities in this codebase).
- Patterns should not change jet geometry (`z_interval`) or bulk flow (`four_velocity`).
"""
abstract type AbstractJetPattern end

"""
    pattern_factor_ne(pattern, x4, jet) -> f

Dimensionless multiplicative factor applied to the *comoving* proper electron density
`electron_density(jet, x4)`.

Return `1` for no modification.

Notes:

- Factors are allowed to be arbitrary (including 0), but the most practical usage
	for physical patterns is typically `1 + …` (a bump/modulation around unity).
"""
function pattern_factor_ne end

"""
    pattern_factor_B(pattern, x4, jet) -> f

Dimensionless multiplicative factor applied to the *comoving* magnetic-field strength
`magnetic_field_strength(jet, x4)`.

Return `1` for no modification.

Notes:

- Factors are allowed to be arbitrary (including 0), but the most practical usage
	for physical patterns is typically `1 + …` (a bump/modulation around unity).
"""
function pattern_factor_B end

@inline pattern_factor_ne(patterns, x4, jet) =
	pattern_factor_ne(only(patterns), x4, jet)
	# prod(p -> pattern_factor_ne(p, x4, jet), patterns)

@inline pattern_factor_B(patterns, x4, jet) =
	prod(p -> pattern_factor_B(p, x4, jet), patterns)

"""Validate that a pattern configuration is physically/semantically allowed for a given jet."""
validate_pattern(::AbstractJetPattern, jet) = nothing


"""
    ConicalBKJetWithPatterns

Wrapper medium that represents:

    (ConicalBKJet baseline) × (multiplicative pattern factors)

Design intent:

- Delegate geometry and kinematics to `base` (same `z_interval`, same `four_velocity`).
- Apply patterns only to comoving `electron_density` and `magnetic_field_strength`.
"""
struct ConicalBKJetWithPatterns{Tbase, Tpatterns} <: AbstractSynchrotronMedium
	base::Tbase
	patterns::Tpatterns
end

function ConicalBKJetWithPatterns(base::ConicalBKJet, patterns)
	for p in patterns
		validate_pattern(p, base)
	end
	return ConicalBKJetWithPatterns{typeof(base), typeof(patterns)}(base, patterns)
end

@unstable prepare_for_computations(obj::ConicalBKJetWithPatterns) = @p let
	obj
	@modify(prepare_for_computations, __.base)
end

@inline z_interval(obj::ConicalBKJetWithPatterns, ray::RayZ) = z_interval(obj.base, ray)
@inline four_velocity(obj::ConicalBKJetWithPatterns, x4) = four_velocity(obj.base, x4)

@inline is_inside_jet(obj::ConicalBKJetWithPatterns, x4::FourPosition) = is_inside_jet(obj.base, x4)

@inline electron_density(obj::ConicalBKJetWithPatterns, x4) =
	electron_density(obj.base, x4) * pattern_factor_ne(obj.patterns, x4, obj.base)

@inline magnetic_field_strength(obj::ConicalBKJetWithPatterns, x4) =
	magnetic_field_strength(obj.base, x4) * pattern_factor_B(obj.patterns, x4, obj.base)

@inline synchrotron_model(obj::ConicalBKJetWithPatterns) = synchrotron_model(obj.base)

abstract type AbstractKnotSizing end

"""
Explicit fixed rest-frame semi-axes `(a_perp, a_parallel)`.

Notes:

- The fields are stored as `(a_perp, a_parallel)`.
- Internally, knot geometry is often handled as `(a_par, a_perp)` to match the
  variable names used in the ellipsoidal form `χ = (Δ_par^2/a_par^2) + (Δ_perp^2/a_perp^2)`.
"""
struct FixedKnotSizing{Ta_perp, Ta_parallel} <: AbstractKnotSizing
	a_perp::Ta_perp
	a_parallel::Ta_parallel
end

"""
Cross-section-tied sizing.

For conical jets the local radius scales as `R(s) = s * tan(φj)`.
This sizing chooses

- `a_perp = f_perp * R(s_c)`
- `a_parallel = q * a_perp`

so the knot occupies a fixed fraction of the local jet cross-section.
"""
struct CrossSectionKnotSizing{Tf_perp, Tq} <: AbstractKnotSizing
	f_perp::Tf_perp
	q::Tq
end

abstract type AbstractKnotProfile end

"""
A smooth bump with `f(χ=0) == f_peak` and `f(χ→∞) → 1`.

Physical meaning:

- If used as `profile_ne`, this can represent a localized overdensity/underdensity
  of radiating electrons comoving with the flow.
- If used as `profile_B`, this can represent a localized magnetic-field enhancement.
"""
@kwdef struct GaussianBump{Tf} <: AbstractKnotProfile
	f_peak::Tf
	χ_threshold::Tf = 4^2
end
GaussianBump(f_peak, χ_threshold) = GaussianBump(promote(f_peak, χ_threshold)...)
GaussianBump(f_peak) = GaussianBump(; f_peak=float(f_peak))

# In the knot frame, χ is an ellipsoidal “radius-squared” and this is a Gaussian in χ.
@inline (g::GaussianBump)(χ) =
	χ < g.χ_threshold ?
		1 + (g.f_peak - 1) * exp(-χ / 2) :
		one(float(g.f_peak))

"""
    InertialEllipsoidalKnot

An SR-consistent “world-tube” pattern with an ellipsoidal profile in the knot rest frame.

- The knot center follows an inertial worldline with constant 4-velocity `u`:
  `x_c(τ) = x_c0 + u * τ`.
- At an event `x4`, we compute the knot proper time at simultaneity in the knot frame
  and the spatial offset in that frame, then evaluate an axisymmetric ellipsoid.

The knot can optionally modulate `n_e` and `|B'|` independently via `profile_ne`/`profile_B`.
"""
@kwdef struct InertialEllipsoidalKnot{TX0, TU, TS, PNE, PB} <: AbstractJetPattern
	x_c0::TX0
	u::TU
	sizing::TS
	profile_ne::PNE
	profile_B::PB
end


_validate_knot_causality(sizing::CrossSectionKnotSizing, u::FourVelocity, jet::ConicalBKJet) = begin
	# Causality guard for the “fills a constant fraction of the cross-section” sizing.
	#
	# With a_perp(τ) = f_perp * R(s_c(τ)) = f_perp * s_c(τ) * tan(φj), we have
	#
	#   v_perp = d a_perp / dτ = f_perp * tan(φj) * d s_c / dτ
	#   v_parallel = d a_parallel / dτ = q * v_perp.
	#
	# For inertial motion, the spatial part of a 4-velocity is u⃗ = d r⃗ / dτ, so
	# d s_c / dτ = jet.axis ⋅ u⃗ (here: ds_dtau = dot(jet.axis, u.xyz)).
	#
	# Physical meaning: if the knot size is tied to the jet radius, the knot must not
	# expand faster than light in its own rest frame.
	ds_dtau = dot(jet.axis, @swiz u.xyz)
	v_perp = abs(sizing.f_perp * tan(jet.φj) * ds_dtau)
	v_parallel = sizing.q * v_perp
	(v_perp < 1 && v_parallel < 1) ||
		error("Unphysical knot sizing: derived rest-frame expansion speeds (v_perp=$v_perp, v_par=$v_parallel) must be < 1")
	return nothing
end

validate_pattern(knot::InertialEllipsoidalKnot, jet::ConicalBKJet) = begin
	knot.sizing isa CrossSectionKnotSizing && _validate_knot_causality(knot.sizing, knot.u, jet)
	return nothing
end


# Return sizes in the order (a_par, a_perp) to match the variables used in `_knot_chi`.
@inline _knot_sizes(sizing::FixedKnotSizing, tau, x_c, jet::ConicalBKJet) = (sizing.a_parallel, sizing.a_perp)

@inline _knot_sizes(sizing::CrossSectionKnotSizing, tau, x_c, jet::ConicalBKJet) = begin
	# Compute rest-frame (a_parallel, a_perp) from the *current* center position.
	# Intuitively: as the knot moves outward, it expands to keep a fixed fraction
	# of the local jet radius.
	s_c = dot(jet.axis, @swiz x_c.xyz)
	a_perp = sizing.f_perp * (s_c * tan(jet.φj))
	a_parallel = sizing.q * a_perp
	return (a_parallel, a_perp)
end


# For a unit timelike u (u⋅u = -1), this is the standard projector onto the local rest space:
#   Δx_perp = (I + u ⊗ u) Δx = Δx + u (u⋅Δx).
# It satisfies u⋅Δx_perp = 0.
@inline _project_orthogonal(u::FourVelocity, dx::FourPosition) = dx + u * minkowski_dot(u, dx)

@inline _knot_velocity_axis(u::FourVelocity) = begin
	# Unit spacelike 4-vector aligned with the knot's lab-frame velocity direction.
	#
	# Math: look for a unit spacelike vector e_par = (a, b β̂) whose spatial direction
	# is along β = beta(u) (β̂ = β/|β|), and which lives in the knot rest space:
	# u⋅e_par = 0 and e_par⋅e_par = 1 (signature (-,+,+,+)).
	β = beta(u)
	βmag = √dot(β, β)
	βhat = β / βmag
	γ = u.t
	return FourPosition(γ * βmag, (γ * βhat)...)
end

@inline _knot_chi(knot::InertialEllipsoidalKnot, x4::FourPosition, jet::ConicalBKJet) = begin
	# Compute the ellipsoidal “radius-squared” χ(x) in the knot rest frame.
	#
	# Overall (metric signature (-,+,+,+)):
	# 1) Find the knot proper time τ such that x_c = x_c0 + u τ is simultaneous with x4
	#    in the knot rest frame (so u⋅(x4 - x_c) = 0).
	# 2) Decompose Δ_rest into components parallel/transverse to the knot principal axis
	#    (taken to be the direction of motion) and evaluate the axisymmetric ellipsoidal
	#    quadratic form.
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

	(a_par, a_perp) = _knot_sizes(knot.sizing, τ, x_c, jet)  # rest-frame semi-axes (parallel, perpendicular)
	return Δ_par^2 / a_par^2 + Δ_perp2 / a_perp^2
end


@inline pattern_factor_ne(knot::InertialEllipsoidalKnot, x4::FourPosition, jet::ConicalBKJet) = begin
	knot.profile_ne === nothing && return one(x4.t)
	χ = _knot_chi(knot, x4, jet)
	return knot.profile_ne(χ)
end

@inline pattern_factor_B(knot::InertialEllipsoidalKnot, x4::FourPosition, jet::ConicalBKJet) = begin
	knot.profile_B === nothing && return one(x4.t)
	χ = _knot_chi(knot, x4, jet)
	return knot.profile_B(χ)
end

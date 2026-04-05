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


# ============================================================================
# GR velocity types for Kerr spacetime
# ============================================================================

# Kerr metric helper functions (M=1 geometrized units, independent of Krang)
_kerr_Δ(a, r) = r^2 - 2r + a^2
_kerr_Σ(a, r, θ) = r^2 + a^2 * cos(θ)^2

# Kerr metric components at (r, θ) in Boyer-Lindquist coordinates
_kerr_g_tt(a, r, θ) = -(1 - 2r / _kerr_Σ(a, r, θ))
_kerr_g_tφ(a, r, θ) = -2a * r * sin(θ)^2 / _kerr_Σ(a, r, θ)
_kerr_g_φφ(a, r, θ) = (r^2 + a^2 + 2a^2 * r * sin(θ)^2 / _kerr_Σ(a, r, θ)) * sin(θ)^2

"""
    KeplerianVelocity(; spin, prograde=true)

Keplerian circular orbit velocity in Kerr spacetime.

At BL radius r (equatorial, M=1):
- Angular velocity: `Ω = ±1 / (r^{3/2} ± a)` (prograde/retrograde)
- 4-velocity: `u^μ = (u^t, 0, 0, Ω·u^t)` with `g_μν u^μ u^ν = -1`
"""
@kwdef struct KeplerianVelocity{T}
	spin::T
	prograde::Bool = true
end

"""
    ZAMOVelocity(; spin)

Zero Angular Momentum Observer in Kerr spacetime.

The ZAMO has `u_φ = 0` (zero covariant angular momentum) and orbits at the
frame-dragging rate `ω = -g_{tφ}/g_{φφ}`. Static at infinity, dragged at finite r.
"""
@kwdef struct ZAMOVelocity{T}
	spin::T
end

"""
    SubKeplerianVelocity(; spin, f_kep, prograde=true)

Fraction of Keplerian angular velocity: `Ω = f_kep · Ω_Keplerian`.
Interpolates between ZAMO-like (f_kep ≈ 0) and Keplerian (f_kep = 1).
"""
@kwdef struct SubKeplerianVelocity{T, Tf}
	spin::T
	f_kep::Tf
	prograde::Bool = true
end


"""
Compute the BL 4-velocity (u^t, 0, 0, u^φ) for a circular orbit with angular velocity Ω
at (r, θ) in Kerr spacetime with spin a.

Returns `(u_t, u_φ)` where `u_t = u^t` and `u_φ = Ω · u^t`.
The normalization `g_μν u^μ u^ν = -1` is enforced.
"""
@inline function _circular_orbit_ut_uφ(a, r, θ, Ω)
	g_tt = _kerr_g_tt(a, r, θ)
	g_tφ = _kerr_g_tφ(a, r, θ)
	g_φφ = _kerr_g_φφ(a, r, θ)
	# g_tt (u^t)² + 2 g_tφ u^t (Ω u^t) + g_φφ (Ω u^t)² = -1
	denom = g_tt + 2g_tφ * Ω + g_φφ * Ω^2
	# denom > 0 means no valid circular orbit (inside ISCO/photon sphere); fall back to rest frame
	-denom ≤ 0 && return (one(Ω), zero(Ω))
	ut = 1 / √(-denom)
	uφ = Ω * ut
	return (ut, uφ)
end

"""
Convert BL 4-velocity (u^t, 0, 0, u^φ) to Cartesian FourVelocity.

The spatial part is pure φ-motion: `v⃗_Cart = (u^φ/u^t) · r sinθ · φ̂_KS`.
The time component `u^t` is preserved as the GR Lorentz factor.
"""
@inline function _bl_circular_to_cartesian_four_velocity(a, r, θ, φ, ut, uφ)
	Ω_coord = uφ / ut  # coordinate angular velocity dφ_BL/dt

	# Convert to KS azimuth for Cartesian direction
	# φ_KS depends on r, so the toroidal direction in KS Cartesian needs the correction.
	# For pure φ-motion (u^r = u^θ = 0), the KS Cartesian velocity is:
	#   v_Cart = Ω_coord · r sinθ · φ̂_KS
	# where φ̂_KS = (-sin φ_KS, cos φ_KS, 0) at the equatorial plane.
	#
	# The ∂φ_KS/∂r correction doesn't matter here because u^r = 0.
	temp = √(1 - a^2)
	rp = 1 + temp
	rm = 1 - temp
	φ_KS = φ - a/(2temp) * log(abs((r - rp)/(r - rm))) + atan(a, r)

	sθ = sin(θ)
	v_φ = Ω_coord * r * sθ  # physical azimuthal speed in coordinate terms
	# φ̂_KS direction in Krang Cartesian:
	φ_hat = SVector(-sin(φ_KS), cos(φ_KS), zero(φ_KS))
	v_cart = v_φ * φ_hat

	# Construct FourVelocity preserving u^t as the GR Lorentz factor
	FourVelocity(ut, ut * v_cart)
end

"""Extract (r, θ, φ) from Cartesian position (assumes BH at origin, spin along z)."""
@inline function _cartesian_to_spherical(xyz::SVector{3})
	r = norm(xyz)
	θ = r > 0 ? acos(clamp(xyz[3] / r, -1, 1)) : zero(r)
	φ = atan(xyz[2], xyz[1])
	return (r, θ, φ)
end

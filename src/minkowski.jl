abstract type FourVector{T} <: FieldVector{4,T} end

struct FourPosition{T} <: FourVector{T}
	t::T
	x::T
	y::T
	z::T
end

struct FourVelocity{T} <: FourVector{T}
    t::T
    x::T
    y::T
    z::T
end

struct FourFrequency{T} <: FourVector{T}
    t::T
    x::T
    y::T
    z::T
end

(::Type{TFV})(t::Number, xyz::SVector{3}) where {TFV<:FourVector} = TFV(t, xyz...)

StaticArrays.similar_type(::Type{TFV}, ::Type{T}, s::Size{(4,)}) where {TFV<:FourVector,T} =
    Base.typename(TFV).wrapper{T}

# lower(v::FourVector) = typeof(v)(-v.t, @swiz v.xyz)
# raise(v::FourVector) = lower(v)

"""
    minkowski_dot(a, b)

Minkowski inner product with signature **(-,+,+,+)**.

For 4-vectors written as `a = (aᵗ, a⃗)` and `b = (bᵗ, b⃗)`:

```
a ⋅ b = −aᵗ bᵗ + a⃗ ⋅ b⃗
```

Used throughout the codebase to avoid manual index raising/lowering.
"""
@inline minkowski_dot(a::FourVector, b::FourVector) = muladd(-a.t, b.t, dot((@swiz a.xyz), (@swiz b.xyz)))

"""
    beta(u)

Return the spatial 3-velocity β⃗ from a timelike 4-vector `u`.

For a 4-velocity written as `u = (γ, γβ⃗)`:

```
β⃗ = u⃗ / uᵗ
```

Errors if `u.t == 0`.
"""
@inline beta(u::FourVelocity) = begin
    iszero(u.t) && error("beta(u): undefined for u.t == 0")
    (@swiz u.xyz) / u.t
end

"""
    direction3(u)

Return the unit spatial direction of motion associated with `u`:

```
v̂ = β⃗ / ‖β⃗‖
```

where `β⃗ = beta(u)`.

This is undefined for a rest 4-velocity (β⃗ = 0).
"""
@inline direction3(u::FourVector) = let
    β = beta(u)
    β / norm(β)
end

@inline direction3(k::FourFrequency) = (@swiz k.xyz) / k.t

"""Unit null 4-vector `(1, n̂) = k/ν`. Used to parameterize ray paths as `x(s) = x₀ + s·direction4(k)`."""
@inline direction4(k::FourFrequency) = k / frequency(k)

"""
    gamma(x)

Lorentz factor γ.

Supported inputs:

- `gamma(u::FourVelocity) = u.t` (expects `uᵗ = γ`).
- `gamma(β::SVector{3}) = 1/√(1 - β⃗⋅β⃗)`.
"""
@inline gamma(u::FourVelocity) = u.t
@inline gamma(β::SVector{3}) = inv(√(1 - dot(β, β)))

# Temporary helper functions for scalar beta <-> gamma conversion
# TODO: Remove once API is standardized to use gamma everywhere
"""
    _beta_from_gamma(γ)

Convert Lorentz factor γ to speed β (temporary helper).
"""
@inline _beta_from_gamma(γ) = √(1 - γ^-2)

"""
    _gamma_from_beta(β)

Convert speed β to Lorentz factor γ (temporary helper).
"""
@inline _gamma_from_beta(β) = inv(√(1 - β^2))

"""
    FourVelocity(β)

Construct a normalized 4-velocity from a spatial 3-velocity β⃗.

```
u = (γ, γβ⃗)
γ = 1 / √(1 - β⃗⋅β⃗)
```
"""
@inline FourVelocity(β::SVector{3}) = let
    γ = gamma(β)
    FourVelocity(γ, γ * β)
end

"""
    construct(FourVelocity, beta=>β, direction3=>dir)

Build a `FourVelocity` from a scalar speed parameter and a spatial unit direction.

- If passed `beta=>β` (a scalar), interprets `β⃗ = β v̂` and returns `FourVelocity(β*dir)`.
- If passed `gamma=>γ` (a scalar), uses `β = √(1 - γ⁻²)` and returns `u = (γ, γβ v̂)`.
"""
@inline construct(::Type{FourVelocity}, (_, β)::Pair{typeof(beta)}, (_, dir)::Pair{typeof(direction3)}) = let
    @boundscheck @assert β isa Number
    @boundscheck @assert dir isa SVector{3}
    return FourVelocity(β * dir)
end

@inline construct(::Type{FourVelocity}, (_, γ)::Pair{typeof(gamma)}, (_, dir)::Pair{typeof(direction3)}) = let
    @boundscheck @assert γ isa Number
    @boundscheck @assert dir isa SVector{3}
	β = _beta_from_gamma(γ)
    return FourVelocity(γ, (γ * β) * dir)
end

"""
    photon_k(ν, n)

Photon 4-frequency in the lab/camera frame for a ray traveling in spatial
direction `n` (unit 3-vector) with frequency ν:

```
k = (ν, ν n)
```
"""
@inline photon_k(νcam, n::SVector{3}) = FourFrequency(νcam, νcam * n)

frequency(k::FourFrequency) = k.t
Accessors.set(k::FourFrequency, ::typeof(frequency), ν) = let
    ν0 = frequency(k)
    return (ν / ν0) * k
end

"""
    comoving_frequency(k, u)

Frequency in the comoving (rest) frame of an observer with 4-velocity `u`:

```
ν′ = −k ⋅ u
```

With the code's metric convention, this is implemented as `-minkowski_dot(k, u)`.
"""
@inline comoving_frequency(k::FourFrequency, u::FourVelocity) = -minkowski_dot(k, u)

"""
    doppler_factor(u, n)

Doppler factor δ for an observer-frame photon direction `n`:

```
δ ≡ ν_obs / ν′
```

In this codebase, `doppler_factor(u,n)` is defined by taking an observer-frame photon
with `ν_obs = 1` and returning

```
δ = 1 / ν′ = 1 / (−k ⋅ u),   with k = (1, n)
```
"""
@inline doppler_factor(u::FourVelocity, n::SVector{3}) = inv(comoving_frequency(photon_k(one(u.t), n), u))

"""
    lorentz_boost_matrix(u)

Return the Lorentz boost matrix `Λ` that transforms a contravariant 4-vector from
the comoving/rest frame of an object into the lab frame in which the object has
4-velocity `u`.

With `u = (γ, γβ⃗)` (so `β⃗ = beta(u)`), this uses the standard boost:

```
Λ⁰₀ = γ
Λ⁰ᵢ = γ βᵢ
Λⁱ₀ = γ βᵢ
Λⁱⱼ = δⁱⱼ + (γ-1) βⁱ βⱼ / (β⃗⋅β⃗)
```

The metric convention is **(-,+,+,+)**.
"""
@inline lorentz_boost_matrix(u::FourVelocity) = let
    γ = u.t
    β = beta(u)
    β2 = dot(β, β)
    iszero(β2) && return @SMatrix [one(γ) zero(γ) zero(γ) zero(γ);
                                   zero(γ) one(γ) zero(γ) zero(γ);
                                   zero(γ) zero(γ) one(γ) zero(γ);
                                   zero(γ) zero(γ) zero(γ) one(γ)]

    A = I + ((γ - 1) / β2) * (β * β')
    γβ = γ * β
    @SMatrix [γ      γβ[1]  γβ[2]  γβ[3];
              γβ[1]  A[1,1] A[1,2] A[1,3];
              γβ[2]  A[2,1] A[2,2] A[2,3];
              γβ[3]  A[3,1] A[3,2] A[3,3]]
end

"""
    lorentz_boost(u, v)

Apply the Lorentz boost defined by `u` to a contravariant 4-vector `v`,
transforming components from the comoving/rest frame into the lab frame.

This is equivalent to `lorentz_boost_matrix(u) * v`, but avoids constructing the
4×4 matrix.
"""
@inline lorentz_boost(u::FourVelocity, v::FourVector) = let
    γ = gamma(u)
    β = beta(u)
    β² = dot(β, β)

    x = @swiz v.xyz
    βx = dot(β, x)

    t′ = muladd(γ, v.t, γ * βx)  # γ*(t + β⋅x)
    s = muladd(γ, v.t, iszero(β) ? 0 : ((γ - 1) / β²) * βx)  # γ*t + ((γ-1)/β²)*(β⋅x)
    x′ = x + s * β
    constructorof(typeof(v))(t′, x′)
end

"""
    lorentz_unboost(u, v)

Apply the inverse Lorentz boost defined by `u` to a contravariant 4-vector `v`,
transforming components from the lab frame (where the object has 4-velocity `u`)
into the object's comoving/rest frame.

This is the inverse of `lorentz_boost(u, ⋅)`.
"""
@inline lorentz_unboost(u::FourVelocity, v::FourVector) = lorentz_boost(FourVelocity(u.t, -@swiz u.xyz), v)

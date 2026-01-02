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

struct FourWavevector{T} <: FourVector{T}
    t::T
    x::T
    y::T
    z::T
end

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
@inline beta(u::FourVector) = begin
    iszero(u.t) && error("beta(u): undefined for u.t == 0")
    (@swiz u.xyz) / u.t
end

"""
    direction(u)

Return the unit spatial direction of motion associated with `u`:

```
v̂ = β⃗ / ‖β⃗‖
```

where `β⃗ = beta(u)`.

This is undefined for a rest 4-velocity (β⃗ = 0).
"""
@inline direction(u::FourVector) = let
    β = beta(u)
    β / norm(β)
end

"""
    gamma(x)

Lorentz factor γ.

Supported inputs:

- `gamma(u::FourVelocity) = u.t` (expects `uᵗ = γ`).
- `gamma(β::SVector{3}) = 1/√(1 - β⃗⋅β⃗)`.
"""
@inline gamma(u::FourVelocity) = u.t
@inline gamma(β::SVector{3}) = inv(√(1 - dot(β, β)))

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
    FourVelocity(γ, (γ * β)...)
end

"""
    construct(FourVelocity, beta=>β, direction=>dir)

Build a `FourVelocity` from a scalar speed parameter and a spatial unit direction.

- If passed `beta=>β` (a scalar), interprets `β⃗ = β v̂` and returns `FourVelocity(β*dir)`.
- If passed `gamma=>γ` (a scalar), uses `β = √(1 - γ⁻²)` and returns `u = (γ, γβ v̂)`.
"""
@inline construct(::Type{FourVelocity}, (_, β)::Pair{typeof(beta)}, (_, dir)::Pair{typeof(direction)}) = let
    @assert β isa Number
    @assert dir isa SVector{3}
    return FourVelocity(β * dir)
end

@inline construct(::Type{FourVelocity}, (_, γ)::Pair{typeof(gamma)}, (_, dir)::Pair{typeof(direction)}) = let
    @assert γ isa Number
    @assert dir isa SVector{3}
	β = √(1 - γ^-2)
    return FourVelocity(γ, (γ * β) * dir...)
end

"""
    photon_k(ν, n)

Photon 4-wavevector (4-frequency) in the lab/camera frame for a ray traveling in spatial
direction `n` (unit 3-vector) with frequency ν:

```
k = (ν, ν n)
```
"""
@inline photon_k(νcam, n::SVector{3}) = FourWavevector(νcam, (νcam .* n)...)

"""
    comoving_frequency(k, u)

Frequency in the comoving (rest) frame of an observer with 4-velocity `u`:

```
ν′ = −k ⋅ u
```

With the code's metric convention, this is implemented as `-minkowski_dot(k, u)`.
"""
@inline comoving_frequency(k::FourWavevector, u::FourVelocity) = -minkowski_dot(k, u)

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
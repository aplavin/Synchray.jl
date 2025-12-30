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

StaticArrays.similar_type(::Type{TFV}, ::Type{T}, s::Size{(4,)}) where {TFV<:FourVector,T} =
    Base.typename(TFV).wrapper{T}

# lower(v::FourVector) = typeof(v)(-v.t, @swiz v.xyz)
# raise(v::FourVector) = lower(v)

@inline minkowski_dot(a::FourVector, b::FourVector) = muladd(-a.t, b.t, dot((@swiz a.xyz), (@swiz b.xyz)))

@inline beta(u::FourVector) = begin
    iszero(u.t) && error("beta(u): undefined for u.t == 0")
    (@swiz u.xyz) / u.t
end
@inline direction(u::FourVector) = let
    β = beta(u)
    β / norm(β)
end

@inline gamma(u::FourVelocity) = u.t
@inline gamma(β::SVector{3}) = inv(√(1 - dot(β, β)))

@inline FourVelocity(β::SVector{3}) = let
    γ = gamma(β)
    FourVelocity(γ, (γ * β)...)
end

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

@inline photon_k(νcam, n::SVector{3}) = FourFrequency(νcam, (νcam .* n)...)

@inline measured_frequency(k::FourFrequency, u::FourVelocity) = -minkowski_dot(k, u)

@inline doppler_factor(u::FourVelocity, n::SVector{3}) = inv(measured_frequency(photon_k(one(u.t), n), u))
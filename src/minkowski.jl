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

lower(v::FourVector) = FourVector(-v.t, @swiz v.xyz)
raise(v::FourVector) = lower(v)

minkowski_dot(a::FourVector, b::FourVector) = -a.t * b.t + dot((@swiz a.xyz), (@swiz b.xyz))

FourVelocity(β::SVector{3}) = let
    γ = 1 / √(1 - dot(β, β))
    FourVelocity(γ, (γ * β)...)
end

photon_k(νcam) = FourFrequency(νcam, zero(νcam), zero(νcam), νcam)

measured_frequency(k::FourFrequency, u::FourVelocity) = -minkowski_dot(k, u)

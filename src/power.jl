import FastPower

struct FixedExponent{P} end
@unstable FixedExponent(p) = FixedExponent{p}()
@inline _half(x::Real) = x / 2
@inline _half(::FixedExponent{x}) where {x} = FixedExponent{x / 2}()

@inline @generated Base.:^(x::Number, ::FixedExponent{P}) where {P} = let
	P == 1.25 && return :(x * √√x)
	isinteger(P) ?
		:(Base.literal_pow(^, x, Val(P))) :
		:(FastPower.fastpower(x, P))
end

Base.:-(a::FixedExponent{P}) where {P} = FixedExponent{-P}()

@inline _value(p::Real) = p
@inline _value(::FixedExponent{P}) where {P} = P



struct StaticNum{X} end

Base.:+(a::Number, ::StaticNum{X}) where {X} = a + X
Base.:+(a::FixedExponent{P}, ::StaticNum{X}) where {P, X} = FixedExponent{P + X}()

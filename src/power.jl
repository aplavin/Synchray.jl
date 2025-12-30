import FastPower

struct FixedExponent{P} end
FixedExponent(p) = FixedExponent{p}()
@inline _half(x::Real) = x / 2
@inline _half(::FixedExponent{x}) where {x} = FixedExponent{x / 2}()

@inline @generated Base.:^(x::Number, ::FixedExponent{P}) where {P} = let
	P == 1.25 && return :(x * √√x)
	isinteger(P) ?
		:(Base.literal_pow(^, x, Val(P))) :
		:(FastPower.fastpower(x, P))
end

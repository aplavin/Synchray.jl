import FastPower

# AD-safe sqrt: forward identical to `sqrt`, but forces the ForwardDiff derivative to 0 at x==0 (its
# true value, since composites like Bperp^q, q>0 vanish there) instead of the spurious inv(2√0)·0 =
# Inf·0 = NaN that poisons the adjoint when B ∥ line of sight (Bperp=0).
@inline sqrt₀(x::Real) = sqrt(x)
@inline function sqrt₀(d::ForwardDiff.Dual{T}) where {T}
	v  = ForwardDiff.value(d)
	sv = sqrt(v)
	deriv = iszero(v) ? zero(sv) : inv(2sv)
	return ForwardDiff.Dual{T}(sv, deriv * ForwardDiff.partials(d))
end

struct FixedExponent{P} end
@unstable FixedExponent(p) = FixedExponent{p}()
@inline _half(x::Real) = x / 2
@inline _half(::FixedExponent{x}) where {x} = FixedExponent{x / 2}()

@inline @generated Base.:^(x::Number, ::FixedExponent{P}) where {P} = let
	P == 1.25 && return :(x * sqrt₀(sqrt₀(x)))
	isinteger(P) && return :(Base.literal_pow(^, x, $(Val(Int(P)))))
	# Half-integer: x^n * √x  (e.g. P=1.5 → x*√x, P=2.5 → x^2*√x, P=0.5 → √x).
	# Uses `sqrt₀` so the AD derivative is finite at x==0 (avoids the √-at-0 `Inf*0=NaN`).
	isinteger(2P) && return :(Base.literal_pow(^, x, $(Val(Int(P - 1//2)))) * sqrt₀(x))
	:(FastPower.fastpower(x, $P))
end

Base.:-(a::FixedExponent{P}) where {P} = FixedExponent{-P}()

@inline _value(p::Real) = p
@inline _value(::FixedExponent{P}) where {P} = P

@unstable to_float_type(::Type{TF}, obj::FixedExponent{P}) where {TF<:AbstractFloat, P} = FixedExponent(TF(P))


struct StaticNum{X} end

Base.:+(a::Number, ::StaticNum{X}) where {X} = a + X
Base.:+(a::FixedExponent{P}, ::StaticNum{X}) where {P, X} = FixedExponent{P + X}()

Base.:*(a::Number, ::StaticNum{X}) where {X} = a * X
Base.:*(a::FixedExponent{P}, ::StaticNum{X}) where {P, X} = FixedExponent{P * X}()

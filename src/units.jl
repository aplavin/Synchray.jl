@kwdef struct UnitsContext{TL,TB,Tn,TF,Tc}
	L0::TL = 1u"cm"        # physical length per 1 code length
	B0::TB = 1u"Gauss"         # physical B-field per 1 code B
	ne0::Tn = 1u"cm^-3"    # physical density per 1 code density
	ν0::TF = 1u"Hz"        # physical frequency per 1 code frequency
	c::Tc = 1u"c"
end

const DEFAULT_CTX = UnitsContext()

_t0(ctx::UnitsContext) = ctx.L0 / ctx.c

_u_to_code(x::Number, scale::AbstractQuantity) = NoUnits(x / scale)
_u_to_code(x::AbstractArray, scale::AbstractQuantity) = _u_to_code.(x, scale)
_u_to_code(x::AbstractInterval, scale::AbstractQuantity) = @modify(x -> _u_to_code(x, scale), endpoints(x)[∗])

withunits(::Type{ConicalBKJet}; ctx::UnitsContext=DEFAULT_CTX, kws...) = let
    kws = NamedTuple(kws)
	ConicalBKJet(;
		kws...,
		φj=NoUnits(kws.φj),
		s=_u_to_code(kws.s, ctx.L0),
		s0=_u_to_code(kws.s0, ctx.L0),
		ne0=_u_to_code(kws.ne0, ctx.ne0),
		B0=_u_to_code(kws.B0, ctx.B0),
	)
end

withunits(::Type{CameraZ}; ctx::UnitsContext=DEFAULT_CTX, kws...) = let
    kws = NamedTuple(kws)
	CameraZ(;
		kws...,
		xys=_u_to_code(kws.xys, ctx.L0),
		ν=_u_to_code(kws.ν, ctx.ν0),
		t=_u_to_code(kws.t, _t0(ctx)),
	)
end


withunits(::typeof(render), cam, obj, what; ctx::UnitsContext=DEFAULT_CTX, kwargs...) =
	_upost_render(render(cam, obj, what; kwargs...), what; ctx)

_upost_render(result, ::Intensity; ctx) = result * 1u"erg/s/cm^2/Hz/sr"

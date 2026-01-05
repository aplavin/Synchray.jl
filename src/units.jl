@kwdef struct UnitsContext{TL,TB,Tn,TF,Tc,TI}
	L0::TL = 1u"cm"        # physical length per 1 code length
	B0::TB = 1u"Gauss"         # physical B-field per 1 code B
	ne0::Tn = 1u"cm^-3"    # physical density per 1 code density
	ν0::TF = 1u"Hz"        # physical frequency per 1 code frequency
	Iν0::TI = 1u"erg/s/cm^2/Hz/sr"  # physical specific intensity per 1 code specific intensity
	c::Tc = 1u"c"
end

const UCTX = UnitsContext()

_t0(ctx::UnitsContext) = ctx.L0 / ctx.c

_u_to_code(x::Number, scale::AbstractQuantity) = NoUnits(x / scale)
_u_to_code(x::AbstractArray, scale::AbstractQuantity) = _u_to_code.(x, scale)
_u_to_code(x::AbstractInterval, scale::AbstractQuantity) = @modify(x -> _u_to_code(x, scale), endpoints(x)[∗])

withunits(::Type{ConicalBKJet}; kws...) = let
    kws = NamedTuple(kws)
	ConicalBKJet(;
		kws...,
		φj=NoUnits(kws.φj),
		s=_u_to_code(kws.s, UCTX.L0),
		s0=_u_to_code(kws.s0, UCTX.L0),
		ne0=_u_to_code(kws.ne0, UCTX.ne0),
		B0=_u_to_code(kws.B0, UCTX.B0),
	)
end

_withunits(pl::PowerLawS, y0_scale) = PowerLawS(pl.exp; val0=_u_to_code(pl.val0, y0_scale), s0=_u_to_code(pl.s0, UCTX.L0))
_withunits(B::BFieldSpec) = BFieldSpec(_withunits(B.scale, UCTX.B0), B.dir, B.wrap)
withunits(::Type{ConicalJet}; kws...) = let
	kws = NamedTuple(kws)
	ConicalJet(;
		kws...,
		φj=NoUnits(kws.φj),
		s=_u_to_code(kws.s, UCTX.L0),
		ne=_withunits(kws.ne, UCTX.ne0),
		B=_withunits(kws.B),
	)
end

withunits(::Type{InertialEllipsoidalKnot}; kws...) = let
	kws = NamedTuple(kws)
	InertialEllipsoidalKnot(;
		kws...,
		x_c0=_u_to_code(kws.x_c0, UCTX.L0),
		# u=_u_to_code(kws.u, UCTX.c),
	)
end

withunits(::Type{CameraZ}; kws...) = let
    kws = NamedTuple(kws)
	CameraZ(;
		kws...,
		xys=_u_to_code(kws.xys, UCTX.L0),
		ν=_u_to_code(kws.ν, UCTX.ν0),
		t=_u_to_code(kws.t, _t0(UCTX)),
	)
end


@unstable withunits(::typeof(render), cam, obj, what; kwargs...) = _upost_render(render(cam, obj, what; kwargs...), what)

_upost_render(result, ::Intensity) = result * UCTX.Iν0
_upost_render(result, ::OpticalDepth) = result  # dimensionless
_upost_render(result, ::SpectralIndex) = result  # dimensionless
_upost_render(result, what::Tuple) = ntuple(length(what)) do i
	_upost_render(getindex.(result, i), what[i])
end

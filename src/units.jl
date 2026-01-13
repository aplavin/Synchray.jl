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


@unstable withunits(::typeof(render), cam, obj, what; kwargs...) = _upost_render(render(cam, obj, what; kwargs...), what)

_upost_render(result, ::Intensity) = result * UCTX.Iν0
_upost_render(result, ::IntensityIQU) = result * UCTX.Iν0
_upost_render(result, ::OpticalDepth) = result  # dimensionless
_upost_render(result, ::SpectralIndex) = result  # dimensionless
_upost_render(result, what::Tuple) = ntuple(length(what)) do i
	_upost_render(getindex.(result, i), what[i])
end

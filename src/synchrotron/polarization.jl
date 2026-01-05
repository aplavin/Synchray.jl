struct StokesIQU{T} <: FieldVector{3,T}
	I::T
	Q::T
	U::T
end

struct ModePerpPar{T} <: FieldVector{2,T}
	perp::T
	par::T
end

struct ModePerpParU{T} <: FieldVector{3,T}
	perp::T
	par::T
	U::T
end

@inline _Pi_j(p) = (p + 1) / (p + 7 / 3)
@inline _Pi_a(p) = (p + 2) / (p + 10 / 3)

@inline modes_from_IQ(I, Q) = ModePerpPar(_half(I + Q), _half(I - Q))
@inline stokes_IQ(m::ModePerpPar) = (I=m.perp + m.par, Q=m.perp - m.par)

@inline emissivity_absorption_polarized(obj::AbstractSynchrotronMedium, x4, k′) = begin
	(jI, αI) = emissivity_absorption(obj, x4, k′)
	field = magnetic_field(obj, x4)

	if field isa FullyTangled
		# Conservative placeholder: unresolved tangling depolarizes.
		j = ModePerpPar(_half(jI), _half(jI))
		a = ModePerpPar(_half(αI), _half(αI))
		return (j, a)
	end
    @assert field isa SVector{3}

	model = synchrotron_model(obj)
	p = _value(model.p)
	Πj = _Pi_j(p)
	Πα = _Pi_a(p)

	j = ModePerpPar((1 + Πj)/2 * jI, (1 - Πj)/2 * jI)
	a = ModePerpPar((1 + Πα)/2 * αI, (1 - Πα)/2 * αI)
	return (j, a)
end

struct PowerLawS{Texp,Tval,Ts0}
	exp::Texp
	val0::Tval
	s0::Ts0
end

PowerLawS(exp; val0, s0=one(val0)) = PowerLawS(exp, val0, s0)

@inline (pl::PowerLawS)(s) = s > 0 ? pl.val0 * (s / pl.s0)^pl.exp : zero(float(pl.val0))

abstract type AbstractJetFieldDirection end

"""Marker direction model meaning: no ordered direction; `BFieldSpec.wrap` receives a scalar strength."""
struct ScalarField <: AbstractJetFieldDirection end

struct PoloidalField <: AbstractJetFieldDirection end
struct ToroidalField <: AbstractJetFieldDirection end
struct HelicalField{T} <: AbstractJetFieldDirection
	ψ::T
end

struct BFieldSpec{Tscale,Tdir,Twrap}
	scale::Tscale
	dir::Tdir
	wrap::Twrap
end

@inline jet_field_direction(::ScalarField, jet, x4) = 1

@inline jet_field_direction(::PoloidalField, jet, x4) = normalize(jet.axis)

@inline jet_field_direction(::ToroidalField, jet, x4) = begin
	axis = normalize(jet.axis)
	r = @swiz x4.xyz
	s = dot(axis, r)
	r_perp = r - s * axis
	ρ = norm(r_perp)
	iszero(ρ) && return zero(axis)
	e_r = r_perp / ρ
	return cross(axis, e_r)
end

@inline jet_field_direction(h::HelicalField, jet, x4) = begin
	axis = normalize(jet.axis)
	r = @swiz x4.xyz
	s = dot(axis, r)
	r_perp = r - s * axis
	ρ = norm(r_perp)
	e_phi = iszero(ρ) ? zero(axis) : cross(axis, r_perp / ρ)
	sψ, cψ = sincos(h.ψ)
	v = cψ * axis + sψ * e_phi
	return iszero(v) ? zero(v) : v / norm(v)
end

@kwdef struct ConicalJet{Ta,Tφ,Ts,Tne,TB,Tu,Tmodel} <: AbstractSynchrotronMedium
	axis::Ta
	φj::Tφ
	s::Ts
	ne::Tne
	B::TB
	speed_profile::Tu
	model::Tmodel
end

@unstable prepare_for_computations(obj::ConicalJet) = @p let
	obj
	@modify(AngleTrigCached_fromangle, __.φj)
	@modify(prepare_for_computations, __.model)
end

z_interval(obj::ConicalJet, ray::RayZ) = _rayz_cone_z_interval(obj.axis, obj.φj, ray, obj.s)

@inline is_inside_jet(jet::ConicalJet, x4::FourPosition) = begin
	r = @swiz x4.xyz
	s = dot(jet.axis, r)
	(s ∈ jet.s) || return false
	ρ = norm(r - s * jet.axis)
	return ρ ≤ s * tan(jet.φj)
end

jet_rotation_matrix(jet::ConicalJet) = jet_rotation_matrix(jet.axis)
jet_basis(jet::ConicalJet) = jet_basis(jet.axis)
lab_to_jet_coords(jet::ConicalJet, r::SVector{3}) = lab_to_jet_coords(jet.axis, r)
jet_to_lab_coords(jet::ConicalJet, rjet::SVector{3}) = jet_to_lab_coords(jet.axis, rjet)

@inline four_velocity(obj::ConicalJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)

	vhat = iszero(r) ? zero(r) : r / √(dot(r, r))
	ρ = norm(r - s * obj.axis)
	tanφ = tan(obj.φj)
	η = ρ / (s * tanφ)

	(whichspeed, speed) = obj.speed_profile(η)
	return construct(FourVelocity, whichspeed => speed, direction => vhat)
end

@inline electron_density(obj::ConicalJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	return obj.ne(s)
end

@inline magnetic_field(obj::ConicalJet, x4) = begin
	r = @swiz x4.xyz
	s = dot(obj.axis, r)
	b = obj.B.scale(s)
	b̂ = jet_field_direction(obj.B.dir, obj, x4)
	return obj.B.wrap(b * b̂)
end

@inline synchrotron_model(obj::ConicalJet) = obj.model

"""
    ConicalJetWithPatterns

Wrapper medium that represents:

    (ConicalJet baseline) × (multiplicative pattern factors)

Design intent:

- Delegate geometry and kinematics to `base`.
- Apply patterns only to comoving `electron_density` and `magnetic_field`.
"""
struct ConicalJetWithPatterns{Tbase, Tpatterns} <: AbstractSynchrotronMedium
	base::Tbase
	patterns::Tpatterns
end

function ConicalJetWithPatterns(base::ConicalJet, patterns)
	for p in patterns
		validate_pattern(p, base)
	end
	return ConicalJetWithPatterns{typeof(base), typeof(patterns)}(base, patterns)
end

@unstable prepare_for_computations(obj::ConicalJetWithPatterns) = @p let
	obj
	@modify(prepare_for_computations, __.base)
end

@inline z_interval(obj::ConicalJetWithPatterns, ray::RayZ) = z_interval(obj.base, ray)
@inline four_velocity(obj::ConicalJetWithPatterns, x4) = four_velocity(obj.base, x4)
@inline is_inside_jet(obj::ConicalJetWithPatterns, x4::FourPosition) = is_inside_jet(obj.base, x4)
@inline jet_basis(obj::ConicalJetWithPatterns) = jet_basis(obj.base)
@inline lab_to_jet_coords(obj::ConicalJetWithPatterns, r::SVector{3}) = lab_to_jet_coords(obj.base, r)
@inline jet_to_lab_coords(obj::ConicalJetWithPatterns, rjet::SVector{3}) = jet_to_lab_coords(obj.base, rjet)

@inline electron_density(obj::ConicalJetWithPatterns, x4) =
	electron_density(obj.base, x4) * pattern_factor_ne(obj.patterns, x4, obj.base)

@inline magnetic_field(obj::ConicalJetWithPatterns, x4) = begin
	field = magnetic_field(obj.base, x4)
	f = pattern_factor_B(obj.patterns, x4, obj.base)
	return field * f
end

@inline synchrotron_model(obj::ConicalJetWithPatterns) = synchrotron_model(obj.base)

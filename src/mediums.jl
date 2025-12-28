abstract type AbstractMedium end

abstract type AbstractSynchrotronMedium <: AbstractMedium end

four_velocity(obj::AbstractMedium, x4, ν) = four_velocity(obj, x4)

# Primary API (what the raytracer calls): invariants
# Default fallbacks are derived from comoving j_ν and α_ν.
emissivity_invariant(obj::AbstractMedium, x4, ν) = emissivity(obj, x4, ν) / (ν^2)
absorption_invariant(obj::AbstractMedium, x4, ν) = absorption(obj, x4, ν) * ν


# --- Synchrotron (Stage 1: angle-averaged Stokes I) ---
@kwdef struct PowerLawElectrons{Tp,TEmin,TEmax,TCj,TCa}
	p::Tp
	Emin::TEmin = nothing
	Emax::TEmax = nothing
	Cj::TCj = 1
	Ca::TCa = 1
end

_synchrotron_coeffs(model::PowerLawElectrons, n_e, B, ν) = begin
	(n_e <= 0 || B <= 0 || ν <= 0) && return (zero(ν), zero(ν))
	(;p, Cj, Ca) = model
	j = Cj * n_e * (B^((p + 1) / 2)) * (ν^(-(p - 1) / 2))
	α = Ca * n_e * (B^((p + 2) / 2)) * (ν^(-(p + 4) / 2))
	return j, α
end

emissivity(obj::AbstractSynchrotronMedium, x4, ν) = begin
	(j, _) = _synchrotron_coeffs(synchrotron_model(obj), electron_density(obj, x4), magnetic_field_strength(obj, x4), ν)
	return j
end

absorption(obj::AbstractSynchrotronMedium, x4, ν) = begin
	(_, α) = _synchrotron_coeffs(synchrotron_model(obj), electron_density(obj, x4), magnetic_field_strength(obj, x4), ν)
	return α
end

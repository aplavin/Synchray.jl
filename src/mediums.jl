abstract type AbstractMedium end

abstract type AbstractSynchrotronMedium <: AbstractMedium end

# Primary API (what the raytracer should prefer calling): compute both comoving j_ν and α_ν
# (they typically share expensive intermediate computations).
emissivity_absorption_invariant(obj::AbstractMedium, x4, ν) = begin
	(j, α) = emissivity_absorption(obj, x4, ν)
	return (j / (ν^2), α * ν)
end

# if a medium defines only `emissivity_absorption`, the generic `emissivity`/`absorption` wrappers below will work
# `absorption` is used in optical depth calculation, maybe can drop this in future?..
emissivity_invariant(obj::AbstractMedium, x4, ν) = emissivity(obj, x4, ν) / (ν^2)
absorption_invariant(obj::AbstractMedium, x4, ν) = absorption(obj, x4, ν) * ν
emissivity(obj::AbstractMedium, x4, ν) = emissivity_absorption(obj, x4, ν)[1]
absorption(obj::AbstractMedium, x4, ν) = emissivity_absorption(obj, x4, ν)[2]


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

emissivity_absorption(obj::AbstractSynchrotronMedium, x4, ν) =
	_synchrotron_coeffs(
		synchrotron_model(obj),
		electron_density(obj, x4), magnetic_field_strength(obj, x4), ν)

abstract type AbstractMedium end

four_velocity(obj::AbstractMedium, x4, ν) = four_velocity(obj, x4)

# Primary API (what the raytracer calls): invariants
# Default fallbacks are derived from comoving j_ν and α_ν.
emissivity_invariant(obj::AbstractMedium, x4, ν) = emissivity(obj, x4, ν) / (ν^2)
absorption_invariant(obj::AbstractMedium, x4, ν) = absorption(obj, x4, ν) * ν

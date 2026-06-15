"""
    BrokenPowerLaw

Broken power-law electron distribution wrapper.

Wraps a single inner electron model (isotropic, anisotropic, or pitchy) for the
low-frequency segment (p₁), and stores only the spectral index p₂ for the
high-frequency segment. The spectrum is continuous at the break by construction:
below ν_break the inner model is evaluated directly, above ν_break the emissivity
and absorption are scaled from their junction values with the p₂ spectral slopes.

The electron distribution is:

    N(γ) = K × (γ/γ_b)^(-p₁)    γ_min ≤ γ < γ_b
         = K × (γ/γ_b)^(-p₂)    γ_b ≤ γ ≤ γ_max

with K = n_e / (I₁ + I₂) determined from the total electron density n_e
(provided by the `ne` profile in `SynchrotronEmission`).

The break frequency is derived from γ_b and the local magnetic field:
`ν_break = (3e)/(4πm_e c) × B × γ_b²`.

# Usage
```julia
# Isotropic broken power law
BrokenPowerLaw(IsotropicPowerLawElectrons; p_low=2.0, p_high=3.0,
    γ_break=Profiles.Transverse(η -> 1e3 + 9e3*η))

# Anisotropic broken power law
BrokenPowerLaw(AnisotropicPowerLawElectrons; p_low=2.0, p_high=3.0, η=0.01,
    γ_break=1e4)  # constant break
```
"""
struct BrokenPowerLaw{Tinner, Tp, Tγb, Tγ}
	inner::Tinner
	p_high::Tp
	γ_break::Tγb
	γmin::Tγ
	γmax::Tγ
end

function BrokenPowerLaw(inner_type; p_low, p_high, γ_break, γmin=1, γmax=Inf, kwargs...)
	@assert p_low > 1 && p_high > 1
	inner = inner_type(; p=p_low, γmin=1, γmax=Inf, kwargs...)
	BrokenPowerLaw(inner, p_high, γ_break, promote(γmin, γmax)...)
end

@unstable function prepare_for_computations(model::BrokenPowerLaw)
	inner = prepare_for_computations(model.inner)
	γb = prepare_for_computations(model.γ_break)
	BrokenPowerLaw(inner, FixedExponent(model.p_high), γb, model.γmin, model.γmax)
end
ustrip(model::BrokenPowerLaw) = model

# Magnetic field strength extraction
@inline _field_strength(f::FullyTangled) = f.strength
@inline _field_strength(b::SVector{3}) = norm(b)
@inline _field_strength(f::TangledOrderedMixture) = norm(f.b)

# K_per_ne for broken power law: n_e = K × (I₁ + I₂), returns 1/(I₁+I₂)
@inline function _K_per_ne_broken(γ_b, γmin, γmax, p_low_raw, p_high_raw)
	p1 = _value(p_low_raw)
	p2 = _value(p_high_raw)
	I1 = γ_b / (p1 - 1) * ((γmin / γ_b)^(1 - p1) - 1)
	I2 = isfinite(γmax) ?
		γ_b / (p2 - 1) * (1 - (γmax / γ_b)^(1 - p2)) :
		γ_b / (p2 - 1)
	return inv(I1 + I2)
end

# Synchrotron break frequency: ν_c(γ_b) = prefactor × B × γ_b²
const _NU_BREAK_PREFACTOR = let
	e  = 4.8032068e-10
	me = ustrip(u"g", u"me")
	c  = ustrip(u"cm/s", 1.0u"c")
	3e / (4π * me * c)
end

@inline function _synchrotron_coeffs(model::BrokenPowerLaw, n_e, field, k′::FourFrequency, geom, x4)
	γ_b = model.γ_break(geom, x4)
	ν = frequency(k′)
	B = _field_strength(field)
	ν_break = oftype(ν, _NU_BREAK_PREFACTOR) * B * γ_b^2

	K_broken = _K_per_ne_broken(γ_b, model.γmin, model.γmax, model.inner.p, model.p_high)
	# Inner model was constructed with γmin=1, γmax=∞ → its K_per_ne = p-1, baked into Cj_ordered.
	# Adjust n_e so that inner's (j,α) = c5(p)*(p-1)*n_e_adj*... gives the correct c5(p)*K_eff*n_e*...
	n_e_adj = n_e * K_broken * γ_b^model.inner.p / (_value(model.inner.p) - 1)

	if ν < ν_break
		return _synchrotron_coeffs(model.inner, n_e_adj, field, k′)
	else
		# Evaluate the inner (low) model at ν_break to get the continuous junction value,
		# then scale with p_high spectral index above the break.
		# j ∝ ν^(-(p-1)/2),  α ∝ ν^(-(p+4)/2)  for power-law synchrotron.
		k_break = k′ * (ν_break / ν)
		(j_b, α_b) = _synchrotron_coeffs(model.inner, n_e_adj, field, k_break)
		ratio = ν / ν_break
		p_hi = model.p_high
		j = j_b * ratio ^ -_half(p_hi + StaticNum{-1}())
		α = α_b * ratio ^ -_half(p_hi + StaticNum{4}())
		return (j, α)
	end
end

"""
		IsotropicPowerLawElectrons(; p, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)

Power-law synchrotron electron model for Stokes-I with **isotropic electron momenta**.

Angle conventions (important):

- theta_Bn: angle between the magnetic-field direction and the line of sight (photon direction)
	in the **plasma rest frame**. This sets `B_perp = |B'| * sin(theta_Bn)`.
- alpha_eB: electron pitch angle, i.e. angle between an electron's velocity and `B'`.

This model does not track electron anisotropy (alpha_eB); the standard power-law synchrotron
coefficients assume an isotropic electron momentum distribution.

Field-direction handling is controlled by the magnetic-field representation passed to
`_synchrotron_coeffs`:

- `FullyTangled(|B′|)` applies an isotropic viewing-angle average over θ_{Bn}.
- An ordered field vector (`SVector{3}`) uses the instantaneous viewing angle.
- `TangledOrderedMixture(b; kappa)` interpolates between these (heuristically).

Normalization:

- If `Cj`/`Ca` are not provided, they are computed from standard cgs synchrotron
	coefficients (Rybicki–Lightman-style `c₅(p)`/`c₆(p)`) combined with the chosen
	power-law normalization, so the model is *physically normalized*.
- If `Cj` and `Ca` are provided explicitly, the model acts as a pure scaling law.

Physical-unit interpretation:

- The core transfer code operates in code units (plain numbers). For physically
	meaningful results with the default normalization, inputs must be consistent with
	cgs conventions for the microphysics (e.g. `n_e` in cm⁻³, `|B′|` in Gauss, `ν` in Hz),
	typically via the Unitful boundary helpers (`withunits`).
"""
struct IsotropicPowerLawElectrons{Tp,Tγ,TC,Tavg}
	p::Tp
	γmin::Tγ
	γmax::Tγ
	Cj_ordered::TC
	Ca_ordered::TC
	sinavg_j::Tavg
	sinavg_a::Tavg
end

@unstable prepare_for_computations(model::IsotropicPowerLawElectrons) = @modify(FixedExponent, model.p)
ustrip(model::IsotropicPowerLawElectrons) = model

# Pitch-angle factor for isotropic electrons: no anisotropy modulation.
@inline _phi_theta(model::IsotropicPowerLawElectrons, cos2θ) = 1

"""
		IsotropicPowerLawElectrons(; p, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)

Constructor for `IsotropicPowerLawElectrons` with **ordered-field** normalization semantics:

- If `Cj`/`Ca` are provided, they are interpreted as the ordered-field coefficients (in terms of `B_perp`).
- If omitted, cgs coefficients are used.
"""
function IsotropicPowerLawElectrons(; p, γmin=1, γmax=Inf, Cj=nothing, Ca=nothing)
	qj = _half(p + 1)
	qa = _half(p + 2)
	sinavg_j = _avg_sin_pow(qj)
	sinavg_a = _avg_sin_pow(qa)

	if isnothing(Cj) || isnothing(Ca)
		@assert isnothing(Cj) && isnothing(Ca)
		K_per_ne = _K_per_ne(p, γmin, γmax)
		(c5, c6) = _synchrotron_c5_c6_ordered(p)
		Cj_ordered = c5 * K_per_ne
		Ca_ordered = c6 * K_per_ne
	else
		# XXX: assume provided Cj/Ca are for tangled field, for backwards compatibility
		Cj_ordered = Cj / sinavg_j
		Ca_ordered = Ca / sinavg_a
	end
	return IsotropicPowerLawElectrons(p, promote(γmin, γmax)..., promote(Cj_ordered, Ca_ordered)..., promote(sinavg_j, sinavg_a)...)
end

# Average of sin^q θ for isotropically distributed directions.
#
# In this codebase's Stage-1 synchrotron, θ is the *viewing angle* θ_{Bn} between the
# magnetic-field direction and the (comoving) photon propagation direction, so
#   B_perp = |B'| sinθ_{Bn}.
#
# If directions are isotropic on the sphere, the relative-angle PDF is p(θ) = (1/2)sinθ on [0, π].
#   ⟨sin^q θ⟩ = (1/2)∫₀^π sin^{q+1}θ dθ
#            = √π · Γ((q+2)/2) / (2 · Γ((q+3)/2)).
@inline _avg_sin_pow(q) = sqrt(pi) * SpecialFunctions.gamma((q + 2) / 2) / (2 * SpecialFunctions.gamma((q + 3) / 2))

# Normalization for a power-law electron distribution N(γ) = K γ^{-p} on [γmin, γmax].
# The number density is
#   n_e = ∫ N(γ) dγ = K/(p-1) · (γmin^(1-p) - γmax^(1-p))  (for p>1).
# Therefore K/n_e = (p-1)/(γmin^(1-p) - γmax^(1-p)); for γmax=Inf, the γmax term → 0.
@inline _K_per_ne(p, γmin, γmax) = begin
	@assert p > 1 && γmin > 0 && γmax > γmin
	invnorm = isfinite(γmax) ? (γmin^(1 - p) - γmax^(1 - p)) : γmin^(1 - p)
	return (p - 1) / invnorm
end

# Analytic coefficient helpers (CGS); see e.g. Rybicki & Lightman.
# Returns the pair (c5(p), c6(p)) used for power-law synchrotron emissivity/absorption.
@inline _synchrotron_c5_c6(p) = begin
	e = 4.8032068e-10           # statC
	me = ustrip(u"g", u"me")     # g
	c = ustrip(u"cm/s", 1.0u"c") # cm/s

	# dimensionless conversion factor from ν/B
	A = (2 * pi * me * c) / (3 * e)
	pref0 = sqrt(3) * e^3 / (16 * pi * me)

	# Stage-1 assumes an isotropically tangled field and averages the viewing-angle factor
	# ⟨sin^q θ_{Bn}⟩ needed for B_perp^q = (|B'| sinθ_{Bn})^q.
	sinavg5 = _avg_sin_pow((p + 1) / 2)
	sinavg6 = _avg_sin_pow((p + 2) / 2)

	c5 = (pref0 / c^2) * (p - 1) * SpecialFunctions.gamma((3p - 1) / 12) * SpecialFunctions.gamma((3p + 7) / 12) * A^(-(p - 1) / 2) * sinavg5
	c6 = pref0 * (p + 2) * SpecialFunctions.gamma((3p + 2) / 12) * SpecialFunctions.gamma((3p + 10) / 12) * A^(-(p + 2) / 2) * sinavg6
	return c5, c6
end

# Same as `_synchrotron_c5_c6`, but for an ordered field without the viewing-angle average.
# In this case, the standard power-law coefficients depend on B_perp explicitly and there is
# no built-in ⟨sin^q θ⟩ factor.
@inline _synchrotron_c5_c6_ordered(p) = begin
	e = 4.8032068e-10           # statC
	me = ustrip(u"g", u"me")     # g
	c = ustrip(u"cm/s", 1.0u"c") # cm/s

	A = (2 * pi * me * c) / (3 * e)
	pref0 = sqrt(3) * e^3 / (16 * pi * me)

	c5 = (pref0 / c^2) * (p - 1) * SpecialFunctions.gamma((3p - 1) / 12) * SpecialFunctions.gamma((3p + 7) / 12) * A^(-(p - 1) / 2)
	c6 = pref0 * (p + 2) * SpecialFunctions.gamma((3p + 2) / 12) * SpecialFunctions.gamma((3p + 10) / 12) * A^(-(p + 2) / 2)
	return c5, c6
end

@inline _synchrotron_coeffs(model::IsotropicPowerLawElectrons, n_e, field::FullyTangled, k′::FourFrequency) = let
	# Stage 1 (angle-averaged) power-law synchrotron, in the comoving frame.
	# Returns (j_ν, α_ν) with ν measured in the plasma rest frame.
	#
	# Implemented scaling:
	#   j_ν = Cj · n_e · B^((p+1)/2) · ν^(-(p-1)/2)
	#   α_ν = Ca · n_e · B^((p+2)/2) · ν^(-(p+4)/2)
	#
	# Notes on normalization:
	# - If the coefficients were auto-derived (see `IsotropicPowerLawElectrons(...)`), then
	#   interpreting `n_e`, `B`, `ν` as (cm⁻³, Gauss, Hz) yields cgs-normalized coefficients.
	# - If `Cj`/`Ca` were provided explicitly, this is a unitless scaling law.
	#
	# Algebraic factorization used below:
	#   common ≡ (B/ν)^(p/2)
	#   j_ν = Cj · n_e · common · √B · √ν
	#   α_ν = Ca · n_e · common · B · ν^-2
	ν = frequency(k′)
	(;p, Cj_ordered, Ca_ordered, sinavg_j, sinavg_a) = model
	Cj_tangled = Cj_ordered * sinavg_j
	Ca_tangled = Ca_ordered * sinavg_a
	B = field.strength
	invν = inv(ν)
	B_over_ν = B * invν
	common = B_over_ν^_half(p)
	j = Cj_tangled * n_e * common * sqrt(B*ν)
	α = Ca_tangled * n_e * common * B * invν^2
	return j, α
end

@inline _synchrotron_coeffs(model::IsotropicPowerLawElectrons, n_e, b::SVector{3}, k′::FourFrequency) = let
	ν = frequency(k′)
	invν = inv(ν)
	(;p, Cj_ordered, Ca_ordered) = model

	# Photon direction in the comoving frame: k′ = (ν, ν n̂) for a null vector.
	n̂ = (@swiz k′.xyz) * invν
	Bperp = norm(cross(b, n̂))

	B_over_ν = Bperp * invν
	common = B_over_ν^_half(p)
	j = Cj_ordered * n_e * common * sqrt(Bperp * ν)
	α = Ca_ordered * n_e * common * Bperp * invν^2
	return j, α
end

@inline _synchrotron_coeffs(model::IsotropicPowerLawElectrons, n_e, field::TangledOrderedMixture, k′::FourFrequency) = let
	ν = frequency(k′)
	(;p, Cj_ordered, Ca_ordered, sinavg_j, sinavg_a) = model
	κ = field.kappa
	@assert κ ≥ 0

	b = field.b
	B = norm(b)
	invν = inv(ν)

	# Ordered viewing angle from the preferred direction.
	n = (@swiz k′.xyz) * invν
	@assert dot(n, n) ≈ 1
	sinθ = norm(cross(b, n)) / B
	sinθ = clamp(sinθ, 0, 1)

	# Minimal ordering model: mix between isotropic-direction average (κ=0) and fully ordered (κ→∞).
	f = κ == Inf ? one(float(κ)) : κ / (one(κ) + κ)

	qj = _half(p + StaticNum{1}())
	qa = _half(p + StaticNum{2}())
	Aj = muladd(f, sinθ^qj - sinavg_j, sinavg_j)
	Aa = muladd(f, sinθ^qa - sinavg_a, sinavg_a)

	B_over_ν = B * invν
	common = B_over_ν^_half(p)
	j = Cj_ordered * n_e * common * sqrt(B * ν) * Aj
	α = Ca_ordered * n_e * common * B * invν^2 * Aa
	return j, α
end

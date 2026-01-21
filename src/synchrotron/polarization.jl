# Linear polarization Stokes vector (I,Q,U) in some chosen screen basis.
#
# Notes:
# - `I` is total specific intensity (or its invariant form, depending on caller).
# - `Q,U` encode linear polarization relative to the basis; changing basis rotates (Q,U) by 2χ.
struct StokesIQU{T} <: FieldVector{3,T}
	I::T
	Q::T
	U::T
end

StaticArrays.similar_type(::Type{<:StokesIQU}, ::Type{T}, s::Size{(3,)}) where {T} = StokesIQU{T}

evpa(s::StokesIQU) = 0.5 * atan(s.U, s.Q)  # electric vector position angle


# Intrinsic synchrotron normal-mode pair (⊥, ∥) in the field-aligned basis.
#
# Convention (Rybicki–Lightman): +Q ≡ I_⊥ − I_∥, where ⊥/∥ refer to the
# E-vector orientation relative to the projected magnetic field B′_⊥.
struct ModePerpPar{T} <: FieldVector{2,T}
	perp::T
	par::T
end


# Mode-space state used by the polarized transfer integrator: (I_⊥, I_∥, U).
#
# In normal-mode space, ⊥ and ∥ decouple for absorption/emission (no Faraday terms).
struct ModePerpParU{T} <: FieldVector{3,T}
	perp::T
	par::T
	U::T
end

# Fractional linear polarization for power-law electrons, for emissivity and absorption:
# Π_j(p) ≡ j_Q / j_I,  Π_α(p) ≡ α_Q / α_I.
@inline _Pi_j(p) = (p + 1) / (p + 7 / 3)
@inline _Pi_a(p) = (p + 2) / (p + 10 / 3)

# Convert between Stokes (I,Q) and intrinsic modes (I_⊥, I_∥):
# I = I_⊥ + I_∥,  Q = I_⊥ − I_∥.
@inline modes_from_IQ(I, Q) = ModePerpPar((I + Q)/2, (I - Q)/2)
@inline modes_from_IQU(s::StokesIQU) = ModePerpParU((s.I + s.Q)/2, (s.I - s.Q)/2, s.U)
@inline stokes_IQ(m::ModePerpPar) = (I=m.perp + m.par, Q=m.perp - m.par)
@inline stokes_IQU(m::ModePerpParU) = StokesIQU(m.perp + m.par, m.perp - m.par, m.U)


# Dispatch on field type for polarized coefficients

@inline _emissivity_absorption_polarized_field(model, jI, αI, field::FullyTangled, k′) = begin
	# Conservative placeholder: unresolved tangling depolarizes.
	j = ModePerpPar(jI/2, jI/2)
	a = ModePerpPar(αI, αI)
	return (j, a)
end

@inline _emissivity_absorption_polarized_field(model::IsotropicPowerLawElectrons, jI, αI, field::TangledOrderedMixture, k′) = begin
	# Incoherent mixture: tangled (unpolarized) + ordered (polarized) components.
	# Effective polarization weighted by intensity contributions.
	# Effective: Π_eff = Π * f * (I_ordered / I_total)

	# Get intrinsic polarization fractions
	p = _value(model.p)
	Πj = _Pi_j(p)
	Πα = _Pi_a(p)

	# Compute mixing fraction
	κ = field.kappa
	f = κ == Inf ? one(float(κ)) : κ / (one(κ) + κ)

	# Recompute viewing angle factors
	ν = frequency(k′)
	b = field.b
	B = norm(b)
	n = (@swiz k′.xyz) / ν
	sinθ = norm(cross(b, n)) / B
	sinθ = clamp(sinθ, 0, 1)

	# Compute angle-factor exponents
	qj = _half(model.p + StaticNum{1}())
	qa = _half(model.p + StaticNum{2}())

	# Compute ordered angle factors
	Aj_ordered = sinθ^qj
	Aa_ordered = sinθ^qa

	# Compute mixed angle factors (as in Stokes-I code)
	Aj_mixed = muladd(f, Aj_ordered - model.sinavg_j, model.sinavg_j)
	Aa_mixed = muladd(f, Aa_ordered - model.sinavg_a, model.sinavg_a)

	# Effective polarization fractions with safe division
	Πj_eff = Aj_mixed > 0 ? Πj * f * Aj_ordered / Aj_mixed : zero(Πj)
	Πα_eff = Aa_mixed > 0 ? Πα * f * Aa_ordered / Aa_mixed : zero(Πα)

	# Split into normal modes using effective fractions
	j = ModePerpPar((1 + Πj_eff)/2 * jI, (1 - Πj_eff)/2 * jI)
	a = ModePerpPar((1 + Πα_eff) * αI, (1 - Πα_eff) * αI)

	return (j, a)
end

@inline _emissivity_absorption_polarized_field(model, jI, αI, field::SVector{3}, k′) = begin
	# Ordered field: use intrinsic polarization fractions.
	p = _value(model.p)
	Πj = _Pi_j(p)
	Πα = _Pi_a(p)

	# This Π-based mode split is exact for the currently implemented synchrotron models.
	# In particular, for `AnisotropicPowerLawElectrons` the direction-dependent factor φ(θBn)
	# multiplies both Stokes-I emissivity and absorption returned by `emissivity_absorption`.
	# Since Π_j = jQ/jI and Π_α = αQ/αI are ratios, this common prefactor cancels and the
	# intrinsic polarization fractions (and thus the mode split) are unchanged.

	j = ModePerpPar((1 + Πj)/2 * jI, (1 - Πj)/2 * jI)
	# For absorption in normal modes: αI = (α_⊥ + α_∥)/2.
	a = ModePerpPar((1 + Πα) * αI, (1 - Πα) * αI)
	return (j, a)
end

@inline emissivity_absorption_polarized(obj::AbstractSynchrotronMedium, x4, k′::FourFrequency) = begin
	# Return *comoving-frame* polarized emissivity/absorption in the intrinsic field basis:
	# (j_⊥, j_∥), (α_⊥, α_∥). For the currently implemented models we split the existing
	# scalar Stokes-I coefficients using the RL thin-limit polarization fractions Π_j, Π_α.
	(jI, αI) = emissivity_absorption(obj, x4, k′)
	field = magnetic_field(obj, x4)
	model = synchrotron_model(obj)
	return _emissivity_absorption_polarized_field(model, jI, αI, field, k′)
end


# Project v into the plane orthogonal to n (i.e. remove the component along n).
@inline _project_to_plane(v::SVector{3}, n::SVector{3}) = v - dot(v, n) * n

@inline _arbitrary_screen_basis(n::SVector{3}) = begin
	# Pick a reference axis not too parallel to n and build a right-handed basis.
	ref = abs(n.x) < 0.9 ? SVector(1, 0, 0) : SVector(0, 1, 0)
	e1 = normalize(cross(n, ref))
	e2 = cross(n, e1)
	return (e1, e2)
end

"""
	comoving_screen_basis(u, k)

Build an orthonormal spatial basis `(e1′, e2′)` spanning the plane orthogonal to the
comoving photon direction `n′`, derived from the fixed camera basis vectors.

Algorithm:

- Unboost the camera basis vectors with `lorentz_unboost(u, ⋅)`.
- Project them into the plane orthogonal to the comoving ray direction.
- Orthonormalize and enforce right-handedness with `e2′ = n′ × e1′`.

Returns `(n′, e1′, e2′)` where all are unit 3-vectors in the comoving frame.
"""
@inline comoving_screen_basis(u::FourVelocity, k::FourFrequency) = begin
	k′ = lorentz_unboost(u, k)
	# Comoving photon direction (unit 3-vector): n′ = k′⃗ / k′ᵗ.
	n′ = direction(k′)
	@assert dot(n′, n′) ≈ 1.0
	e1_cam = SVector(1, 0, 0)
	e2_cam = SVector(0, 1, 0)

	# Unboost camera basis vectors (treated as spatial four-vectors),
	# then project them into the screen plane ⟂ n′.
	e1′raw = @swiz lorentz_unboost(u, FourPosition(0, e1_cam)).xyz
	e2′raw = @swiz lorentz_unboost(u, FourPosition(0, e2_cam)).xyz

	e1′p = _project_to_plane(e1′raw, n′)
	e1′ = e1′p / norm(e1′p)
	e2′ = cross(n′, e1′)

	# Align the sign with the projected unboosted camera e2 (a 180° flip is Stokes-invariant).
	e2′p = _project_to_plane(e2′raw, n′)
	if dot(e2′, e2′p) < 0
		e1′ = -e1′
		e2′ = -e2′
	end
	@assert !iszero(e1′p) && !iszero(e2′p)

	return (n′, e1′, e2′)
end

"""
	linear_polarization_basis_from_B(n′, B′)

Construct a right-handed comoving screen basis aligned to the projected magnetic field.

Returns `(e_par, e_perp)` where:

- `e_par` is the unit vector along the projection of `B′` into the plane ⟂ `n′`.
- `e_perp = n′ × e_par`.

Convention note (used throughout the polarized transfer code):

- The *field-aligned Stokes basis* takes the +Q axis along `e_perp` (electric field
	perpendicular to the projected magnetic field).
- With this choice, `ModePerpPar(perp, par)` corresponds to intensities with the
	E-vector along `(e_perp, e_par)` respectively, and in the field basis
	Q ≡ I_⊥ − I_∥.

If `B′ ∥ n′` (vanishing projection), uses an arbitrary screen basis instead.
"""
@inline linear_polarization_basis_from_B(n′::SVector{3}, B′::SVector{3}) = begin
	# Return both screen-plane axes tied to B′: e_par ∥ proj(B′), e_perp ⟂ proj(B′).
	b̂ = iszero(B′) ? B′ : normalize(B′)
	@assert all(!isnan, b̂)
	bp = _project_to_plane(b̂, n′)
	if dot(bp, bp) == 0
		(e_par, e_perp) = _arbitrary_screen_basis(n′)
		return (e_par, e_perp)
	end
	e_par = bp / norm(bp)
	e_perp = cross(n′, e_par)
	return (e_par, e_perp)
end

linear_polarization_basis_from_B(n′::SVector{3}, B′::FullyTangled) = (SVector(1,0,0), SVector(0,1,0))

linear_polarization_basis_from_B(n′::SVector{3}, B′::TangledOrderedMixture) = linear_polarization_basis_from_B(n′, B′.b)

@inline stokes_QU_rotation(χ) = let
	# Rotate Stokes (Q,U) for a basis rotation by χ in the screen plane (2χ law).
	s, c = sincos(χ)
	c2 = muladd(c, c, -(s * s))
	s2 = 2 * s * c
	@SMatrix [c2 s2; -s2 c2]
end

"""
	stokes_QU_rotation(e1_old, e2_old, e1_new)

Rotation matrix acting on the Stokes `(Q,U)` components when changing the linear
polarization basis from `(e1_old, e2_old)` to `(e1_new, e2_new)` in the *same* screen plane.

The returned 2×2 matrix `R` satisfies:

```
[Q; U]_new = R * [Q; U]_old
```

and uses the standard 2χ law for Stokes parameters.
	
"""
@inline stokes_QU_rotation(e1_old::SVector{3}, e2_old::SVector{3}, e1_new::SVector{3}) = begin
	# Compute cosχ,sinχ via dot products (no explicit angle extraction), then apply 2χ law.
	cosχ = dot(e1_new, e1_old)
	sinχ = dot(e1_new, e2_old)
	c2 = muladd(cosχ, cosχ, -(sinχ * sinχ))
	s2 = 2 * sinχ * cosχ
	@SMatrix [c2 s2; -s2 c2]
end

# Apply a Q/U rotation matrix to a (Q,U) 2-vector.
@inline rotate_QU(R, QU::SVector{2}) = R * QU

# Rotate full Stokes (I,Q,U); I is invariant under basis rotations.
@inline rotate_QU(R, s::StokesIQU) = let
	QU = rotate_QU(R, @swiz s.QU)
	StokesIQU(s.I, QU...)
end

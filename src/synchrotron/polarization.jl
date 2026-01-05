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

@inline modes_from_IQ(I, Q) = ModePerpPar((I + Q)/2, (I - Q)/2)
@inline stokes_IQ(m::ModePerpPar) = (I=m.perp + m.par, Q=m.perp - m.par)

@inline emissivity_absorption_polarized(obj::AbstractSynchrotronMedium, x4, k′) = begin
	(jI, αI) = emissivity_absorption(obj, x4, k′)
	field = magnetic_field(obj, x4)

	if field isa FullyTangled
		# Conservative placeholder: unresolved tangling depolarizes.
		j = ModePerpPar(jI/2, jI/2)
		a = ModePerpPar(αI/2, αI/2)
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
	n′ = (@swiz k′.xyz) / k′.t
	@assert dot(n′, n′) ≈ 1.0
	e1_cam = SVector(1, 0, 0)
	e2_cam = SVector(0, 1, 0)

	e1′raw = @swiz lorentz_unboost(u, FourPosition(k′.t, e1_cam)).xyz
	e2′raw = @swiz lorentz_unboost(u, FourPosition(k′.t, e2_cam)).xyz

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

If `B′ ∥ n′` (vanishing projection), uses an arbitrary screen basis instead.
"""
@inline linear_polarization_basis_from_B(n′::SVector{3}, B′::SVector{3}) = begin
	b̂ = normalize(B′)
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

@inline stokes_QU_rotation(χ) = let
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
	cosχ = dot(e1_new, e1_old)
	sinχ = dot(e1_new, e2_old)
	c2 = muladd(cosχ, cosχ, -(sinχ * sinχ))
	s2 = 2 * sinχ * cosχ
	@SMatrix [c2 s2; -s2 c2]
end

@inline rotate_QU(R, QU::SVector{2}) = R * QU

@inline rotate_IQU(R, s::StokesIQU) = let
	(Q, U) = rotate_QU(R, s.Q, s.U)
	StokesIQU(s.I, Q, U)
end

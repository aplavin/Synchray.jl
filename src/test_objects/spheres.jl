"""General ray-sphere intersection: returns interval in ray parameter `s`."""
_sphere_z_interval(center_xyz, radius, ray::Ray) = begin
	n̂ = ray_direction(ray)
	r0 = @swiz ray.x0.xyz
	Δ = center_xyz - r0
	s_c = dot(Δ, n̂)
	b2 = dot(Δ, Δ) - s_c^2
	r2 = radius^2
	ds = b2 > r2 ? -eps(float(s_c)) : √(r2 - b2)
	return (s_c - ds) .. (s_c + ds)
end

@kwdef struct UniformSphere{TC,TR,TU,FJ,FA} <: AbstractMedium
	center::TC
	radius::TR
	u0::TU
	jν::FJ   # comoving emissivity j_ν (constant)
	αν::FA   # comoving absorption α_ν (constant)
end

z_interval(obj::UniformSphere, ray::Ray) =
	_sphere_z_interval(@swiz(obj.center.xyz), obj.radius, ray)

four_velocity(obj::UniformSphere, x4) = obj.u0

emissivity_absorption(obj::UniformSphere, x4, k′) = (obj.jν, obj.αν)



@kwdef struct MovingUniformEllipsoid{TC,TSZ,TU,FJ,FA} <: AbstractMedium
	center::TC
	sizes::TSZ
	u0::TU
	jν::FJ   # comoving emissivity j_ν (constant)
	αν::FA   # comoving absorption α_ν (constant)
end

_spatial_in_rest(u, v) = begin
	β = beta(u)
	β2 = dot(β, β)
	xyz = (@swiz v.xyz)
	if iszero(β2)
		return xyz
	end
	γ = u.t
	f = ((γ - 1) * dot(β, xyz) / β2 - γ * v.t)
	return xyz + f * β
end

z_interval(obj::MovingUniformEllipsoid, ray::Ray) = begin
	# Worldtube of a rigid axis-aligned ellipsoid moving with constant 4-velocity u.
	# Ray: x(s) = ray.x0 + s·k1, where k1 = k/ν = (1, n̂).
	u = obj.u0
	Δ0 = ray.x0 - obj.center
	ν = frequency(ray)
	kdir = ray.k / ν  # = (1, n̂)

	a = minkowski_dot(u, kdir)
	b = minkowski_dot(u, Δ0)
	P0 = Δ0 + u * b
	P1 = kdir + u * a

	p0 = _spatial_in_rest(u, P0)
	p1 = _spatial_in_rest(u, P1)

	s² = obj.sizes .^ 2
	A = sum(p1 .^ 2 ./ s²)
	B = 2 * sum(p0 .* p1 ./ s²)
	C = sum(p0 .^ 2 ./ s²) - one(A)
	D = B^2 - 4 * A * C

	if !(D > 0) || iszero(A)
		s0 = zero(A)
		return s0 .. (s0 - eps(oneunit(s0)))
	end

	sD = √(D)
	s1 = (-B - sD) / (2 * A)
	s2 = (-B + sD) / (2 * A)
	return min(s1, s2) .. max(s1, s2)
end

four_velocity(obj::MovingUniformEllipsoid, x4) = obj.u0
emissivity_absorption(obj::MovingUniformEllipsoid, x4, k′) = (obj.jν, obj.αν)


@kwdef struct UniformSynchrotronSphere{TC,TR,TU,Tne,TB,TM} <: AbstractMedium
	center::TC
	radius::TR
	u0::TU
	ne0::Tne
	B0::TB
	electrons::TM
end

z_interval(obj::UniformSynchrotronSphere, ray::Ray) =
	_sphere_z_interval(@swiz(obj.center.xyz), obj.radius, ray)

four_velocity(obj::UniformSynchrotronSphere, x4) = obj.u0

electron_density(obj::UniformSynchrotronSphere, x4) = obj.ne0
magnetic_field(obj::UniformSynchrotronSphere, x4) = FullyTangled(obj.B0)
emissivity_absorption(obj::UniformSynchrotronSphere, x4, k′) =
	_synchrotron_coeffs(obj.electrons, obj.ne0, FullyTangled(obj.B0), k′)
emissivity_absorption_polarized(obj::UniformSynchrotronSphere, x4, k′) = begin
	(jI, αI) = emissivity_absorption(obj, x4, k′)
	(j, α) = _emissivity_absorption_polarized_field(obj.electrons, jI, αI, FullyTangled(obj.B0), k′)
	return (j, α, FullyTangled(obj.B0))
end

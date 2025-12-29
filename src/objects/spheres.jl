@kwdef struct UniformSphere{TC,TR,TU,FJ,FA} <: AbstractMedium
	center::TC
	radius::TR
	u0::TU
	jν::FJ   # comoving emissivity j_ν (constant)
	αν::FA   # comoving absorption α_ν (constant)
end

z_interval(obj::UniformSphere, ray::RayZ) = begin
	z0 = obj.center.z
	r2 = obj.radius^2
	dxy = (@swiz ray.x0.xy) - (@swiz obj.center.xy)
	b2 = dot(dxy, dxy)
	dz = b2 > r2 ? -eps(float(z0)) : √(r2 - b2)
	return (z0 - dz) .. (z0 + dz)
end

four_velocity(obj::UniformSphere, x4) = obj.u0

emissivity(obj::UniformSphere, x4, ν) = obj.jν
absorption(obj::UniformSphere, x4, ν) = obj.αν



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

z_interval(obj::MovingUniformEllipsoid, ray::RayZ) = begin
	# Worldtube of a rigid axis-aligned ellipsoid moving with constant 4-velocity u.
	# In the comoving frame, define Δ⊥ as the displacement orthogonal to u, and require
	# (Δx/sx)^2 + (Δy/sy)^2 + (Δz/sz)^2 ≤ 1.
	# `sizes == SVector(R,R,R)` recovers the moving sphere.
	# For RayZ: x(z) = ray.x0 + z*(1,0,0,1) (since k/kz = (1,0,0,1)).
	u = obj.u0
	Δ0 = ray.x0 - obj.center
	onez = one(Δ0.t)
	zeroz = zero(Δ0.t)
	kdir = FourPosition(onez, zeroz, zeroz, onez)

	a = minkowski_dot(u, kdir)
	b = minkowski_dot(u, Δ0)
	P0 = Δ0 + u * b
	P1 = kdir + u * a

	p0 = _spatial_in_rest(u, P0)
	p1 = _spatial_in_rest(u, P1)

	invsizes2 = inv.(obj.sizes .^ 2)
	A = sum(invsizes2 .* (p1 .^ 2))
	B = 2 * sum(invsizes2 .* (p0 .* p1))
	C = sum(invsizes2 .* (p0 .^ 2)) - one(A)
	D = B^2 - 4 * A * C

	if !(D > 0) || iszero(A)
		z0 = obj.center.z
		return z0 .. (z0 - eps(float(z0)))
	end

	√D = sqrt(D)
	z1 = (-B - √D) / (2 * A)
	z2 = (-B + √D) / (2 * A)
	return min(z1, z2) .. max(z1, z2)
end

four_velocity(obj::MovingUniformEllipsoid, x4) = obj.u0
emissivity(obj::MovingUniformEllipsoid, x4, ν) = obj.jν
absorption(obj::MovingUniformEllipsoid, x4, ν) = obj.αν


@kwdef struct UniformSynchrotronSphere{TC,TR,TU,Tne,TB,TM} <: AbstractSynchrotronMedium
	center::TC
	radius::TR
	u0::TU
	ne0::Tne
	B0::TB
	model::TM
end

z_interval(obj::UniformSynchrotronSphere, ray::RayZ) = begin
	z0 = obj.center.z
	r2 = obj.radius^2
	dxy = (@swiz ray.x0.xy) - (@swiz obj.center.xy)
	b2 = dot(dxy, dxy)
	dz = b2 > r2 ? -eps(float(z0)) : √(r2 - b2)
	return (z0 - dz) .. (z0 + dz)
end

four_velocity(obj::UniformSynchrotronSphere, x4) = obj.u0

electron_density(obj::UniformSynchrotronSphere, x4) = obj.ne0
magnetic_field_strength(obj::UniformSynchrotronSphere, x4) = obj.B0
synchrotron_model(obj::UniformSynchrotronSphere) = obj.model

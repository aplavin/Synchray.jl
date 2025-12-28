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



@kwdef struct MovingUniformSphere{TC,TR,TU,FJ,FA} <: AbstractMedium
	center::TC
	radius::TR
	u0::TU
	jν::FJ   # comoving emissivity j_ν (constant)
	αν::FA   # comoving absorption α_ν (constant)
end

z_interval(obj::MovingUniformSphere, ray::RayZ) = begin
	# Worldtube of a rigid sphere moving with constant 4-velocity u:
	# for an event x, define Δ = x - center and project it onto the hyperplane
	# orthogonal to u:  Δ⊥ = Δ + u (u⋅Δ). Points satisfy Δ⊥⋅Δ⊥ ≤ R^2.
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

	A = minkowski_dot(P1, P1)
	B = 2 * minkowski_dot(P0, P1)
	C = minkowski_dot(P0, P0) - obj.radius^2
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

four_velocity(obj::MovingUniformSphere, x4) = obj.u0
emissivity(obj::MovingUniformSphere, x4, ν) = obj.jν
absorption(obj::MovingUniformSphere, x4, ν) = obj.αν


@kwdef struct UniformSlab{TZ,TU,TJ,TA} <: AbstractMedium
	z::TZ
	u0::TU
	j0::TJ
	a0::TA
end

z_interval(obj::UniformSlab, ray::RayZ) = obj.z
four_velocity(obj::UniformSlab, x4) = obj.u0
emissivity(obj::UniformSlab, x4, ν) = obj.j0
absorption(obj::UniformSlab, x4, ν) = obj.a0


@kwdef struct UniformSynchrotronSlab{TZ,TU,Tne,TB,TM} <: AbstractSynchrotronMedium
	z::TZ
	u0::TU
	ne0::Tne
	B0::TB
	model::TM
end

z_interval(obj::UniformSynchrotronSlab, ray::RayZ) = obj.z
four_velocity(obj::UniformSynchrotronSlab, x4) = obj.u0

electron_density(obj::UniformSynchrotronSlab, x4) = obj.ne0
magnetic_field_strength(obj::UniformSynchrotronSlab, x4) = obj.B0
synchrotron_model(obj::UniformSynchrotronSlab) = obj.model

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

four_velocity(obj::UniformSphere, x4, ν) = obj.u0

four_velocity(obj::UniformSphere, x4) = obj.u0

emissivity(obj::UniformSphere, x4, ν) = obj.jν
absorption(obj::UniformSphere, x4, ν) = obj.αν


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

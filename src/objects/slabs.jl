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

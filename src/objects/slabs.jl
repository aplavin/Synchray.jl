@kwdef struct UniformSlab{TZ,TU,TJ,TA} <: AbstractMedium
	z::TZ
	u0::TU
	j0::TJ
	a0::TA
end

z_interval(obj::UniformSlab, ray::RayZ) = obj.z
four_velocity(obj::UniformSlab, x4) = obj.u0
emissivity_absorption(obj::UniformSlab, x4, k′) = (obj.j0, obj.a0)


@kwdef struct UniformSynchrotronSlab{TZ,TU,Tne,TB,TM} <: AbstractMedium
	z::TZ
	u0::TU
	ne0::Tne
	B0::TB
	electrons::TM
end

z_interval(obj::UniformSynchrotronSlab, ray::RayZ) = obj.z
four_velocity(obj::UniformSynchrotronSlab, x4) = obj.u0

electron_density(obj::UniformSynchrotronSlab, x4) = obj.ne0
magnetic_field(obj::UniformSynchrotronSlab, x4) = obj.B0
emissivity_absorption(obj::UniformSynchrotronSlab, x4, k′) =
	_synchrotron_coeffs(obj.electrons, obj.ne0, obj.B0, k′)
emissivity_absorption_polarized(obj::UniformSynchrotronSlab, x4, k′) = begin
	(jI, αI) = emissivity_absorption(obj, x4, k′)
	(j, α) = _emissivity_absorption_polarized_field(obj.electrons, jI, αI, obj.B0, k′)
	return (j, α, obj.B0)
end

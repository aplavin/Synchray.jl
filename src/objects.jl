@kwdef struct UniformSphere{FJ,FA} <: AbstractMedium
	center::FourPosition
	radius
	u0::FourVelocity
	jν::FJ   # comoving emissivity j_ν(ν) or constant-like callable
	αν::FA   # comoving absorption α_ν(ν) or constant-like callable
end

z_interval(obj::UniformSphere, x, y, νcam, t0) = begin
	x0, y0, z0 = obj.center[2], obj.center[3], obj.center[4]
	dx, dy = x - x0, y - y0
	r2 = obj.radius^2
	b2 = dx^2 + dy^2
	b2 > r2 ? nothing : begin
		dz = sqrt(r2 - b2)
		(z0 - dz) .. (z0 + dz)
	end
end

four_velocity(obj::UniformSphere, x4, ν) = obj.u0

four_velocity(obj::UniformSphere, x4) = obj.u0

emissivity(obj::UniformSphere, x4, ν) = obj.jν(ν)
absorption(obj::UniformSphere, x4, ν) = obj.αν(ν)


@kwdef struct ConstantSlab <: AbstractMedium
	z::Interval
	u0::FourVelocity
	j0
	a0
end

z_interval(obj::ConstantSlab, x, y, νcam, t0) = obj.z
four_velocity(obj::ConstantSlab, x4) = obj.u0
emissivity(obj::ConstantSlab, x4, ν) = obj.j0
absorption(obj::ConstantSlab, x4, ν) = obj.a0

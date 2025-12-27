@kwdef struct UniformSphere{TC,TR,TU,FJ,FA} <: AbstractMedium
	center::TC
	radius::TR
	u0::TU
	jν::FJ   # comoving emissivity j_ν (constant)
	αν::FA   # comoving absorption α_ν (constant)
end

z_interval(obj::UniformSphere, xy::SVector{2}, νcam, t0) = begin
	x0, y0, z0 = obj.center[2], obj.center[3], obj.center[4]
	dx, dy = xy[1] - x0, xy[2] - y0
	r2 = obj.radius^2
	b2 = dx^2 + dy^2
	b2 > r2 ? z0..(z0 - eps(float(z0))) : begin
		dz = √(r2 - b2)
		(z0 - dz) .. (z0 + dz)
	end
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

z_interval(obj::UniformSlab, xy, νcam, t0) = obj.z
four_velocity(obj::UniformSlab, x4) = obj.u0
emissivity(obj::UniformSlab, x4, ν) = obj.j0
absorption(obj::UniformSlab, x4, ν) = obj.a0

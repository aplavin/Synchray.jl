@kwdef struct OrthoCamera
	xys
    nz::Int
    ν
    t
    mapfunc = map
end

@unstable render(cam::OrthoCamera, obj::AbstractMedium) = begin
	map(cam.xys) do xy::SVector{2}
		integrate_ray(obj, xy; νcam=cam.ν, t0=cam.t, cam.nz)
	end
end

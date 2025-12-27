render(cam::OrthoCamera, obj::AbstractMedium; νcam, t0=0.0, nz=256) = begin
	rows = map(cam.ys) do y
		map(cam.xs) do x
			integrate_ray(obj, x, y; νcam, t0, nz)
		end
	end
	reduce(vcat, permutedims.(rows))
end

struct OrthoCamera
	xs
	ys
end

OrthoCamera(xmin, xmax, nx, ymin, ymax, ny) = OrthoCamera(range(xmin, xmax; length=nx), range(ymin, ymax; length=ny))

module AdaptExt

using Synchray: CombinedMedium, CameraOrtho
import Adapt

Adapt.adapt_structure(to, cm::CombinedMedium) = CombinedMedium(Adapt.adapt(to, cm.objects))

Adapt.adapt_structure(to, cam::CameraOrtho) =
	CameraOrtho(
		cam.origin, cam.n, cam.e1, cam.e2,
		Adapt.adapt(to, cam.xys),
		cam.nz, cam.ν, cam.t, cam.light,
		map,
	)

end

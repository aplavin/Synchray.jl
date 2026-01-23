@testitem "Adaptive supersampling: boundary pixels" begin
	import Synchray as S
	using RectiGrids

	sphere = S.UniformSphere(
		center=S.FourPosition(0, 0, 0, 0),
		radius=1,
		u0=S.FourVelocity(SVector(0, 0, 0)),
		jν=1,
		αν=0,
	)

	# Grid chosen so one pixel center lies near the boundary but still inside.
	xs = range(-1.5..1.5, step=0.3)  # step == 0.3, includes x=0.9 and y=0
	cam = S.CameraZ(; xys=grid(SVector, xs, xs), nz=512, ν=1, t=0)

	img0 = S.render(cam, sphere)
	img3 = S.render(cam, sphere; adaptive_supersampling=3)
	@test !(img3 ≈ img0)

	Icenter = img0(0, 0)
	@test img3(0, 0) ≈ Icenter

	# Edge-adjacent pixel should be reduced due to partial coverage.
	Iedge0 = img0(0.9, 0)
	Iedge3 = img3(0.9, 0)
	@test Iedge0 > 0
	@test Iedge3 > 0
	@test Iedge3 / Iedge0 < 0.9

	# Fully outside remains zero.
	@test img0(1.5, 0) == 0
	@test img3(1.5, 0) == 0
end

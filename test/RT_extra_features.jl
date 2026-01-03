
@testitem "Ray contribution profiles" begin
	import Synchray as S

	# Simple, constant-coefficient medium: UniformSlab.
	L = 3.0
	j0 = 0.7
	a0 = 1.3
	slab = S.UniformSlab(0.0..L, S.FourVelocity(SVector(0.0, 0.0, 0.0)), j0, a0)

	ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=2.0, nz=2048)

	prof = S.ray_contribution_profile(slab, ray)
	Iν = S.render(ray, slab, S.Intensity())
	τ = S.render(ray, slab, S.OpticalDepth())

	@test Iν > 0
	@test τ > 0
	@test sum(prof.Δτ) ≈ τ rtol=2e-3
	@test sum(prof.dIν_to_obs) ≈ Iν rtol=2e-3

	@testset "missed rays return empty/zero" begin
		R = 1.3
		sphere = S.UniformSphere(; center=S.FourPosition(0.0, 0.0, 0.0, 0.0), radius=R, u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)), jν=1.0, αν=2.0)
		miss_ray = S.RayZ(; x0=S.FourPosition(0.0, 2R, 0.0, 0.0), k=2.0, nz=256)

		prof_miss = S.ray_contribution_profile(sphere, miss_ray)
		@test isempty(prof_miss.z)
		@test isempty(prof_miss.Δτ)
		@test isempty(prof_miss.dIν_to_obs)

		@test S.render(miss_ray, sphere, S.Intensity()) == 0
		@test S.render(miss_ray, sphere, S.OpticalDepth()) == 0
	end
end

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

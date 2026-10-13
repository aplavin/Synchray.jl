@testitem "I contribution profiles" begin
	import Synchray as S

	# Simple, constant-coefficient medium: UniformSlab.
	L = 3.0
	j0 = 0.7
	a0 = 1.3
	slab = S.UniformSlab(0.0..L, S.FourVelocity(SVector(0.0, 0.0, 0.0)), j0, a0)

	ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=2.0, nz=2048)

	prof = S.ray_contribution_profile(slab, ray)
	Iν = @inferred S.render(ray, slab, S.Intensity())
	τ = S.render(ray, slab, S.OpticalDepth())

	# stationary slab, δ=1: τ = a0·L = 3.9; Iν = (j0/a0)(1-e^{-a0 L})
	@test Iν ≈ 0.528 rtol=1e-2
	@test τ ≈ 3.9 rtol=1e-2
	@test sum(prof.Δτ) ≈ τ rtol=2e-3
	@test sum(prof.dIν_to_obs) ≈ Iν rtol=2e-3

	@testset "Midpoint cells, nz=$nz" for nz in (1, 2, 8)
		ray_cells = @set ray.nz = nz
		profile = S.ray_contribution_profile(slab, ray_cells)
		@test length(profile) == nz
		@test all(s -> 0 < s < L, profile.s)
		@test sum(profile.Δs) ≈ L
		@test sum(profile.Δτ) ≈ a0 * L
		@test first(profile.s) ≈ L / (2nz)
		@test last(profile.s) ≈ L - L / (2nz)
		@test sum(profile.dIν_to_obs) ≈ S.render(ray_cells, slab)
	end

	@testset "missed rays return empty/zero" begin
		R = 1.3
		sphere = S.UniformSphere(; center=S.FourPosition(0.0, 0.0, 0.0, 0.0), radius=R, u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)), jν=1.0, αν=2.0)
		miss_ray = S.RayZ(; x0=S.FourPosition(0.0, 2R, 0.0, 0.0), k=2.0, nz=256)

		prof_miss = S.ray_contribution_profile(sphere, miss_ray)
		@test isempty(prof_miss.s)
		@test isempty(prof_miss.Δτ)
		@test isempty(prof_miss.dIν_to_obs)

		@test S.render(miss_ray, sphere, S.Intensity()) == 0
		@test S.render(miss_ray, sphere, S.OpticalDepth()) == 0
	end
end

@testitem "ray_contribution_profile_IQU uniform slab" begin
	import Synchray as S

	# Baseline configuration: uniform slab with ordered magnetic field
	L = 2.0
	B_ordered = SVector(0.0, 0.5, 0.0)  # Ordered field in y direction
	p = 2.5
	electrons = S.IsotropicPowerLawElectrons(; p=p, γmin=10.0, γmax=1e4, Cj=1e-3, Ca=1e-5)

	make_medium(; z=0..L, u0=S.FourVelocity(SVector(0,0,0)), B0=B_ordered, electrons=electrons) =
		S.UniformSynchrotronSlab(; z, u0, ne0=1.0, B0, electrons)

	medium = make_medium()
	ray = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2.0, nz=1_000)

	@testset "Conservation: sum equals full integration" begin
		S_direct = @inferred S.render(ray, medium, S.IntensityIQU())
		profile = S.ray_contribution_profile_IQU(medium, ray)
		S_sum = sum(profile.dIν_to_obs)

		@test S_sum ≈ S_direct rtol=1e-5
	end

	@testset "Scalar consistency: I matches ray_contribution_profile" begin
		profile_scalar = S.ray_contribution_profile(medium, ray)
		profile_IQU = S.ray_contribution_profile_IQU(medium, ray)
		@test [step.dIν_to_obs.I for step in profile_IQU] ≈ profile_scalar.dIν_to_obs rtol=1e-5
	end

	@testset "Physical validity" begin
		profile = S.ray_contribution_profile_IQU(medium, ray)
		for step in profile
			@test step.dIν_to_obs.I ≥ 0
			pol_frac_sq = (step.dIν_to_obs.Q^2 + step.dIν_to_obs.U^2) / (step.dIν_to_obs.I^2 + 1e-20)
			@test pol_frac_sq ≤ 1.0 + 1e-6
		end
	end

	@testset "Midpoint cells, nz=$nz" for nz in (1, 2, 8)
		ray_min = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2.0, nz)
		profile = S.ray_contribution_profile_IQU(medium, ray_min)
		profile_scalar = S.ray_contribution_profile(medium, ray_min)
		@test length(profile) == nz
		@test profile.s == profile_scalar.s
		@test profile.Δs == profile_scalar.Δs
		@test sum(profile.Δs) ≈ L

		S_direct = S.render(ray_min, medium, S.IntensityIQU())
		S_sum = sum(profile.dIν_to_obs)
		@test S_sum ≈ S_direct rtol=1e-3
	end

	@testset "FullyTangled: zero polarization" begin
		medium_tangled = make_medium(; B0=S.FullyTangled(0.5))
		ray_tangled = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2.0, nz=500)

		S_direct = S.render(ray_tangled, medium_tangled, S.IntensityIQU())
		profile = S.ray_contribution_profile_IQU(medium_tangled, ray_tangled)
		S_sum = sum(profile.dIν_to_obs)

		@test S_direct.I ≈ 3.54e-4 rtol=1e-2
		@test abs(S_direct.Q) < 1e-12
		@test abs(S_direct.U) < 1e-12
		@test S_sum ≈ S_direct rtol=1e-5
	end
end


@testitem "ray_contribution_profile_IQU complex geometry" begin
	import Synchray as S

	# EmissionRegion + Conical geometry + Helical field
	# Field orientation varies along the ray path
	θ = 0.1
	region = S.EmissionRegion(
		geometry = S.Geometries.Conical(; axis=SVector(sin(θ), 0, cos(θ)), φj=0.1, z=0.5..20.0),
		velocity = S.VelocitySpec(S.Directions.Axial(), S.gamma, S.Profiles.Constant(10)),
		emission = S.SynchrotronEmission(
			ne = S.Profiles.Constant(1.0),
			B = S.BFieldSpec(
				S.Profiles.Constant(0.5),
				S.Directions.HelicalAT(π/6),
				identity
			),
			electrons = S.IsotropicPowerLawElectrons(; p=2.5, Cj=2e-3, Ca=5e-5),
		),
	)

	ray = S.RayZ(; x0=S.FourPosition(0, 1, 0, 0), k=2.0, nz=500)

	@testset "Conservation with spatially varying field" begin
		profile = S.ray_contribution_profile_IQU(region, ray)
		S_sum = sum(profile.dIν_to_obs)
		S_direct = S.render(ray, region, S.IntensityIQU())

		@test all(!=(0), S_direct)
		@test S_sum ≈ S_direct rtol=3e-3
	end

	@testset "Empty ray returns empty profile" begin
		# Ray far off-axis that misses the conical geometry
		ray_miss = S.RayZ(; x0=S.FourPosition(0, 100, 0, 0), k=2.0, nz=256)
		profile = S.ray_contribution_profile_IQU(region, ray_miss)
		@test isempty(profile)
	end
end
@testitem "step_contribution" begin
	import Synchray as S

	# stationary ball with constant (j, α) = (0.7, 1.3); ν = 2, Δλ = 0.05 ⇒ Δτ = 0.13 (exponential branch)
	ball = S.EmissionRegion(geometry=S.Geometries.Ball(S.Geometries.StationaryAtOrigin(), 5.0), emission=S.FixedEmission(S=0.7 / 1.3, α=1.3))
	st = S.RayState(S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=2.0, nz=1))
	c = @inferred S.step_contribution(ball, S.FourPosition(0.0, 0.0, 0.0, 2.0), st, 0.05)
	@test c.Δτ ≈ 0.13
	@test c.dIinv_source ≈ 0.06564092181200544
	# equals one rendering step from zero intensity
	@test c.dIinv_source ≈ S._integrate_ray_step(S.AccValue{S.Intensity}(0.0), ball, S.FourPosition(0.0, 0.0, 0.0, 2.0), st, 0.05).value
end

@testitem "attenuate_to_observer" begin
	import Synchray as S

	p = S.attenuate_to_observer(S.StructArray(; s=[1, 2, 3], Δτ=[0.1, 0.2, 0.3], dIinv_source=[1.0, 2.0, 3.0]))
	@test p.s == [1, 2, 3]
	@test p.τ_front ≈ [0.5, 0.3, 0.0]
	@test p.dIν_to_obs ≈ [0.6065306597126334, 1.4816364413634358, 3.0]
end

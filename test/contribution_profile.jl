@testitem "I contribution profiles" begin
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
		S_direct = S.render(ray, medium, S.IntensityIQU())
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

	@testset "Minimal steps (nz=2)" begin
		ray_min = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2.0, nz=2)
		profile = S.ray_contribution_profile_IQU(medium, ray_min)
		@test length(profile) == 2

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

		@test S_direct.I > 0
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

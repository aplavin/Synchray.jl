@testitem "EmissionRegion with Knot patterns" begin
	import Synchray as S
	using RectiGrids
	using Accessors

	# Base emission region (equivalent to the ConicalJet from patterns_knots.jl)
	region = S.EmissionRegion(
		geometry = S.Geometries.Conical(; axis = SVector(0, 0, 1), φj = 0.05, z = 1e-3..10),
		ne = S.Profiles.Axial(S.PowerLaw(-2; val0=2., s0=1.)),
		B = S.BFieldSpec(
			S.Profiles.Axial(S.PowerLaw(-1; val0=3., s0=1.)),
			S.Directions.Scalar(),
			b -> S.FullyTangled(b),
		),
		velocity = S.VelocitySpec(S.Directions.Axial(), S.Profiles.Constant(1.0)),
		model = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1, Ca=1),
	)

	# A moving knot centered on-axis at z=2
	knot = S.Patterns.EllipsoidalKnot(
		x_c0 = S.FourPosition(0.0, 0.0, 0.0, 2.0),
		u = S.FourVelocity(SVector(0.0, 0.0, 0.1)),
		sizing = S.Patterns.FixedSizing(0.2, 0.5),
		profile = S.Patterns.GaussianBump(2.0),
	)

	# Create region with knot patterns
	region_with_knot = @set region.ne = S.Profiles.Modified(region.ne, knot)
	region_with_knot = @set region_with_knot.B.scale = S.Profiles.Modified(region_with_knot.B.scale, knot)

	x4c = knot.x_c0
	@test knot(region.geometry, x4c, 1.0) ≈ 2.0
	@test S.electron_density(region_with_knot, x4c) ≈ 2 * S.electron_density(region, x4c)
	@test S.magnetic_field(region_with_knot, x4c) ≈ 2 * S.magnetic_field(region, x4c)

	@testset "knot center stays peaked at different lab times" begin
		center_event_at_lab_time(t) = begin
			τ = (t - knot.x_c0.t) / knot.u.t
			knot.x_c0 + knot.u * τ
		end

		for t in (0.0, 1.0, 5.0)
			x4ct = center_event_at_lab_time(t)
			@test x4ct.t ≈ t
			@test S._knot_chi(knot, region.geometry, x4ct) ≈ 0
			@test knot(region.geometry, x4ct, 1.0) ≈ 2.0
			@test S.electron_density(region_with_knot, x4ct) ≈ 2 * S.electron_density(region, x4ct)
			@test S.magnetic_field(region_with_knot, x4ct) ≈ 2 * S.magnetic_field(region, x4ct)
		end
	end

	@testset "off-center transverse offset gives intermediate factors" begin
		# For motion along +z and a purely transverse offset Δ = (0, a_perp, 0, 0) at the same lab time,
		# the construction in `_knot_chi` yields Δ_par = 0 and χ = (a_perp^2) / (a_perp^2) = 1.
		t = 5.0
		τ = (t - knot.x_c0.t) / knot.u.t
		x4ct = knot.x_c0 + knot.u * τ
		a_perp = knot.sizing.a_perp
		x4off = x4ct + S.FourPosition(0.0, a_perp, 0.0, 0.0)
		@test S._knot_chi(knot, region.geometry, x4off) ≈ 1 atol=5e-12

		f_ne_expected = 1 + (2.0 - 1) * exp(-1 / 2)
		@test knot(region.geometry, x4off, 1.0) ≈ f_ne_expected
		@test S.electron_density(region_with_knot, x4off) ≈ f_ne_expected * S.electron_density(region, x4off)
		@test S.magnetic_field(region_with_knot, x4off) ≈ f_ne_expected * S.magnetic_field(region, x4off)
	end

	# Far from the knot: factors should approach 1.
	x4far = S.FourPosition(0.0, 0.0, 0.0, 5.0)
	@test knot(region.geometry, x4far, 1.0) ≈ 1.0 rtol=0 atol=5e-8
	@test S.electron_density(region_with_knot, x4far) ≈ S.electron_density(region, x4far) rtol=0 atol=5e-8
	@test S.magnetic_field(region_with_knot, x4far) ≈ S.magnetic_field(region, x4far) rtol=0 atol=5e-8

	# Geometry and flow must remain delegated to the base region.
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=128)
	@test S.z_interval(region_with_knot, ray) == S.z_interval(region, ray)
	@test S.four_velocity(region_with_knot, x4c) == S.four_velocity(region, x4c)

	@testset "float64 vs float32" begin
		cam64 = S.CameraZ(; xys=grid(S.SVector, range(0.01..0.1, 2), range(-0.001..0.001, 2)), nz=20, ν=2., t=0.)
		cam32 = S.to_float_type(Float32, cam64)

		region64 = region_with_knot
		region32 = S.to_float_type(Float32, S.prepare_for_computations(region64))

		f64 = S.render(cam64, region64)
		f32 = S.render(cam32, region32)
		@test eltype(f64) == Float64
		@test eltype(f32) == Float32
		@test all(>(0), f64)
		@test f32 ≈ f64  rtol=1e-5
	end

	@testset "EmissionRegion renders match with prepare_for_computations" begin
		regionp = S.prepare_for_computations(region_with_knot)
		cam = S.CameraZ(; xys=grid(S.SVector, range(0.01..0.1, 6), range(-0.001..0.001, 6)), nz=50, ν=2., t=0.)

		img_region = S.render(cam, region_with_knot)
		img_regionp = S.render(cam, regionp)
		@test !(img_regionp ≈ S.render(cam, region))
		@test img_regionp ≈ img_region
	end

	@testset "CrossSectionSizing knot" begin
		# Test knot with CrossSectionSizing
		knot_cs = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 2.0),
			u = S.FourVelocity(SVector(0.0, 0.0, 0.1)),
			sizing = S.Patterns.CrossSectionSizing(0.3, 2.0),
			profile = S.Patterns.GaussianBump(2.5),
		)

		region_cs = @set region.ne = S.Profiles.Modified(region.ne, knot_cs)
		
		# At the center, should see the peak factor
		x4c = knot_cs.x_c0
		@test S.electron_density(region_cs, x4c) ≈ 2.5 * S.electron_density(region, x4c)
		
		# Check that sizes scale with position
		# At z=2, a_perp should be 0.3 * 2 * tan(0.05) ≈ 0.0300
		z_c = 2.0
		expected_a_perp = 0.3 * z_c * tan(0.05)
		expected_a_par = 2.0 * expected_a_perp
		(a_par, a_perp) = S._knot_sizes(knot_cs.sizing, 0.0, x4c, region.geometry)
		@test a_perp ≈ expected_a_perp
		@test a_par ≈ expected_a_par
	end

	@testset "Causality validation" begin
		# Valid knot should pass
		knot_valid = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 1.0),
			u = construct(S.FourVelocity, S.gamma=>5, S.direction=>SVector(0,0,1)),
			sizing = S.Patterns.CrossSectionSizing(0.1, 2.0),
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test S.validate_pattern(knot_valid, region.geometry) === nothing

		# Invalid knot with superluminal expansion should fail
		knot_invalid = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 1.0),
			u = construct(S.FourVelocity, S.gamma=>10, S.direction=>SVector(0,0,1)),
			sizing = S.Patterns.CrossSectionSizing(0.8, 30.0),  # Large q causes superluminal expansion
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test_throws ErrorException S.validate_pattern(knot_invalid, region.geometry)

		# FixedSizing doesn't need causality validation
		knot_fixed = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 1.0),
			u = construct(S.FourVelocity, S.gamma=>10, S.direction=>SVector(0,0,1)),
			sizing = S.Patterns.FixedSizing(0.1, 0.2),
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test S.validate_pattern(knot_fixed, region.geometry) === nothing
	end
end

@testitem "Jet patterns: wrapper + inertial knot" begin
	import Synchray as S
	using RectiGrids
	using Accessors

	jet = S.ConicalJet(;
		axis=SVector(0., 0, 1),
		φj=0.05,
		s=1e-3..10,
		ne=S.PowerLawS(-2; val0=2., s0=1.),
		B=S.BFieldSpec(S.PowerLawS(-1; val0=3., s0=1.), S.ScalarField(), b -> S.FullyTangled(b)),
		speed_profile=(η -> (S.beta, 0)),
		model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1, Ca=1),
	)

	# A moving knot centered on-axis at s=2.
	knot = S.InertialEllipsoidalKnot(
		x_c0=S.FourPosition(0.0, 0.0, 0.0, 2.0),
		u=S.FourVelocity(SVector(0.0, 0.0, 0.1)),
		sizing=S.FixedKnotSizing(0.2, 0.5),
		profile_ne=S.GaussianBump(2.0),
		profile_B=S.GaussianBump(3.0),
	)

	jetp = S.JetWithPatterns(jet, (knot,))

	x4c = knot.x_c0
	@test S.pattern_factor_ne(knot, x4c, jet) ≈ 2.0
	@test S.pattern_factor_B(knot, x4c, jet) ≈ 3.0
	@test S.electron_density(jetp, x4c) ≈ 2 * S.electron_density(jet, x4c)
	@test S.magnetic_field(jetp, x4c) ≈ 3 * S.magnetic_field(jet, x4c)

	@testset "knot center stays peaked at different lab times" begin
		center_event_at_lab_time(t) = begin
			τ = (t - knot.x_c0.t) / knot.u.t
			knot.x_c0 + knot.u * τ
		end

		for t in (0.0, 1.0, 5.0)
			x4ct = center_event_at_lab_time(t)
			@test x4ct.t ≈ t
			@test S._knot_chi(knot, x4ct, jet) ≈ 0
			@test S.pattern_factor_ne(knot, x4ct, jet) ≈ 2.0
			@test S.pattern_factor_B(knot, x4ct, jet) ≈ 3.0
			@test S.electron_density(jetp, x4ct) ≈ 2 * S.electron_density(jet, x4ct)
			@test S.magnetic_field(jetp, x4ct) ≈ 3 * S.magnetic_field(jet, x4ct)
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
		@test S._knot_chi(knot, x4off, jet) ≈ 1 atol=5e-12

		f_ne_expected = 1 + (2.0 - 1) * exp(-1 / 2)
		f_B_expected = 1 + (3.0 - 1) * exp(-1 / 2)
		@test S.pattern_factor_ne(knot, x4off, jet) ≈ f_ne_expected
		@test S.pattern_factor_B(knot, x4off, jet) ≈ f_B_expected
		@test S.electron_density(jetp, x4off) ≈ f_ne_expected * S.electron_density(jet, x4off)
		@test S.magnetic_field(jetp, x4off) ≈ f_B_expected * S.magnetic_field(jet, x4off)
	end

	# Far from the knot: factors should approach 1.
	x4far = S.FourPosition(0.0, 0.0, 0.0, 5.0)
	@test S.pattern_factor_ne(knot, x4far, jet) ≈ 1.0 rtol=0 atol=5e-8
	@test S.pattern_factor_B(knot, x4far, jet) ≈ 1.0 rtol=0 atol=5e-8
	@test S.electron_density(jetp, x4far) ≈ S.electron_density(jet, x4far) rtol=0 atol=5e-8
	@test S.magnetic_field(jetp, x4far) ≈ S.magnetic_field(jet, x4far) rtol=0 atol=5e-8

	# Geometry and flow must remain delegated to the base jet.
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=128)
	@test S.z_interval(jetp, ray) == S.z_interval(jet, ray)
	@test S.four_velocity(jetp, x4c) == S.four_velocity(jet, x4c)

	@testset "float64 vs float32" begin
		cam64 = S.CameraZ(; xys=grid(S.SVector, range(0.01..0.1, 2), range(-0.001..0.001, 2)), nz=20, ν=2., t=0.)
		cam32 = S.to_float_type(Float32, cam64)

		jet64 = jetp
		jet32 = S.to_float_type(Float32, S.prepare_for_computations(jet64))

		f64 = S.render(cam64, jet64)
		f32 = S.render(cam32, jet32)
		@test eltype(f64) == Float64
		@test eltype(f32) == Float32
		@test all(>(0), f64)
		@test f32 ≈ f64  rtol=1e-5
	end

	@testset "ConicalJet renders match for patterns with prepare_for_computations" begin
		cjp = S.JetWithPatterns(jet, (knot,))
		cjpp = S.prepare_for_computations(cjp)
		cam = S.CameraZ(; xys=grid(S.SVector, range(0.01..0.1, 6), range(-0.001..0.001, 6)), nz=50, ν=2., t=0.)

		img_jet = S.render(cam, jetp)
		img_cjp = S.render(cam, cjp)
		img_cjpp = S.render(cam, cjpp)
		@test !(img_cjp ≈ S.render(cam, jet))
		@test img_cjp ≈ img_jet
		@test img_cjpp ≈ img_jet
	end
end

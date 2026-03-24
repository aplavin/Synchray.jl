@testitem "EmissionRegion with Knot patterns" begin
	import Synchray as S
	using RectiGrids
	using Accessors

	# Base emission region (equivalent to the ConicalJet from patterns_knots.jl)
	region = S.EmissionRegion(
		geometry = S.Geometries.Conical(; axis = SVector(0, 0, 1), φj = 0.05, z = 1e-3..10),
		velocity = S.VelocitySpec(S.Directions.Axial(), S.Profiles.Constant(1.0)),
		emission = S.SynchrotronEmission(
			ne = S.Profiles.Axial(S.PowerLaw(-2; val0=2., s0=1.)),
			B = S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=3., s0=1.)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
			electrons = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1, Ca=1),
		),
	)

	# A moving knot centered on-axis at z=2
	knot = S.Patterns.EllipsoidalKnot(
		x_c0 = S.FourPosition(0.0, 0.0, 0.0, 2.0),
		u = S.FourVelocity(SVector(0.0, 0.0, 0.1)),
		sizing = S.Patterns.FixedSizing(0.2, 0.5),
		profile = S.Patterns.GaussianBump(2.0),
	)

	# Create region with knot patterns
	region_with_knot = @set region.emission.ne = S.Profiles.Modified(region.emission.ne, knot)
	region_with_knot = @set region_with_knot.emission.B.scale = S.Profiles.Modified(region_with_knot.emission.B.scale, knot)

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

		region_cs = @set region.emission.ne = S.Profiles.Modified(region.emission.ne, knot_cs)
		
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
			u = construct(S.FourVelocity, S.gamma=>5, S.direction3=>SVector(0,0,1)),
			sizing = S.Patterns.CrossSectionSizing(0.1, 2.0),
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test S.validate_pattern(knot_valid, region.geometry) === nothing

		# Invalid knot with superluminal expansion should fail
		knot_invalid = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 1.0),
			u = construct(S.FourVelocity, S.gamma=>10, S.direction3=>SVector(0,0,1)),
			sizing = S.Patterns.CrossSectionSizing(0.8, 30.0),  # Large q causes superluminal expansion
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test_throws ErrorException S.validate_pattern(knot_invalid, region.geometry)

		# FixedSizing doesn't need causality validation
		knot_fixed = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 1.0),
			u = construct(S.FourVelocity, S.gamma=>10, S.direction3=>SVector(0,0,1)),
			sizing = S.Patterns.FixedSizing(0.1, 0.2),
			profile = S.Patterns.GaussianBump(2.0),
		)
		@test S.validate_pattern(knot_fixed, region.geometry) === nothing
	end

end

@testitem "retarded_event" begin
	import Synchray as S
	using Accessors
	using RectiGrids

	# Geometry needed for _knot_chi
	geom = S.Geometries.Conical(; axis=SVector(0, 0, 1), φj=0.05, z=1e-3..10)

	# Knot moving along z-axis (simple case)
	knot_axial = S.Patterns.EllipsoidalKnot(
		x_c0 = S.FourPosition(0.0, 0.0, 0.0, 2.0),
		u = S.FourVelocity(SVector(0.0, 0.0, 0.5)),
		sizing = S.Patterns.FixedSizing(0.2, 0.5),
		profile = S.Patterns.GaussianBump(2.0),
	)

	# Knot moving diagonally (off-axis, relativistic)
	knot_diag = S.Patterns.EllipsoidalKnot(
		x_c0 = S.FourPosition(0.0, 1.0, 0.5, 3.0),
		u = construct(S.FourVelocity, S.gamma => 5, S.direction3 => normalize(SVector(0.2, 0.1, 0.97))),
		sizing = S.Patterns.FixedSizing(0.1, 0.1),
		profile = S.Patterns.GaussianBump(2.0),
	)

	# Default Z-camera for tests
	cam_z = S.CameraZ(; xys=grid(SVector, x=[0.0], y=[0.0]), nz=1, ν=1.0, t=0.0)

	@testset "worldline method" begin
		x0 = S.FourPosition(0.0, 1.0, 0.5, 3.0)
		u = S.FourVelocity(SVector(0.1, 0.05, 0.9))
		wl = S.Geometries.InertialWorldline(x0, u)
		for t_obs in (0.0, 1.0, -1.0)
			cam = @set cam_z.t = t_obs
			ret = S.retarded_event(wl, cam)
			# Null hyperplane condition: t - z = t_obs
			@test ret.x.t - ret.x.z ≈ t_obs atol=1e-12
			# Verify on worldline: x = x0 + u*tau
			@test ret.x ≈ x0 + u * ret.tau atol=1e-12
		end
	end

	@testset "knot method matches worldline method" for knot in (knot_axial, knot_diag)
		wl = S.Geometries.InertialWorldline(knot.x_c0, knot.u)
		for t_obs in (0.0, 2.0)
			cam = @set cam_z.t = t_obs
			ret_knot = S.retarded_event(knot, cam)
			ret_wl = S.retarded_event(wl, cam)
			@test ret_knot.x ≈ ret_wl.x atol=1e-12
			@test ret_knot.tau ≈ ret_wl.tau atol=1e-12
		end
	end

	@testset "consistency with _knot_chi" for knot in (knot_axial, knot_diag)
		for t_obs in (0.0, 1.0, 5.0, -2.0)
			cam = @set cam_z.t = t_obs
			ret = S.retarded_event(knot, cam)
			# ret.x is the knot center FourPosition; chi should be 0 there
			@test S._knot_chi(knot, geom, ret.x) ≈ 0 atol=1e-10
		end
	end

	@testset "round-trip with camera_ray_anchor" for knot in (knot_axial, knot_diag)
		for t_obs in (0.0, 1.0, 5.0)
			cam = @set cam_z.t = t_obs
			ret = S.retarded_event(knot, cam)
			# camera_ray_anchor should recover (x, y) and t_obs
			anchor = S.camera_ray_anchor(cam, ret.x)
			@test anchor.x ≈ ret.x.x atol=1e-12
			@test anchor.y ≈ ret.x.y atol=1e-12
			@test anchor.t ≈ t_obs atol=1e-12
		end
	end

	@testset "null hyperplane condition" for knot in (knot_axial, knot_diag)
		for t_obs in (0.0, 1.0, -1.0)
			cam = @set cam_z.t = t_obs
			ret = S.retarded_event(knot, cam)
			# Retarded event lies on null hyperplane t - z = t_obs
			@test ret.x.t - ret.x.z ≈ t_obs atol=1e-12
		end
	end

	@testset "FastLight mode" begin
		cam_fast_z = S.CameraZ(; xys=grid(SVector, x=[0.0], y=[0.0]), nz=1, ν=1.0, t=0.0, light=S.FastLight())

		# FastLight worldline: simultaneity condition t = t_obs
		@testset "worldline simultaneity" begin
			x0 = S.FourPosition(0.0, 1.0, 0.5, 3.0)
			u = S.FourVelocity(SVector(0.1, 0.05, 0.9))
			wl = S.Geometries.InertialWorldline(x0, u)
			for t_obs in (0.0, 1.0, -1.0)
				cam = @set cam_fast_z.t = t_obs
				ret = S.retarded_event(wl, cam)
				# FastLight: event is at t = t_obs (simultaneity, not null hyperplane)
				@test ret.x.t ≈ t_obs atol=1e-12
				# Verify on worldline
				@test ret.x ≈ x0 + u * ret.tau atol=1e-12
			end
		end

		# FastLight knot method matches worldline method
		@testset "knot matches worldline" for knot in (knot_axial, knot_diag)
			wl = S.Geometries.InertialWorldline(knot.x_c0, knot.u)
			for t_obs in (0.0, 2.0)
				cam = @set cam_fast_z.t = t_obs
				ret_knot = S.retarded_event(knot, cam)
				ret_wl = S.retarded_event(wl, cam)
				@test ret_knot.x ≈ ret_wl.x atol=1e-12
				@test ret_knot.tau ≈ ret_wl.tau atol=1e-12
			end
		end

		# FastLight: _knot_chi = 0 at retarded event center
		@testset "knot_chi consistency" for knot in (knot_axial, knot_diag)
			for t_obs in (0.0, 1.0, 5.0, -2.0)
				cam = @set cam_fast_z.t = t_obs
				ret = S.retarded_event(knot, cam)
				@test S._knot_chi(knot, geom, ret.x) ≈ 0 atol=1e-10
			end
		end

		# FastLight: round-trip with camera_ray_anchor
		@testset "camera_ray_anchor roundtrip" for knot in (knot_axial, knot_diag)
			for t_obs in (0.0, 1.0, 5.0)
				cam = @set cam_fast_z.t = t_obs
				ret = S.retarded_event(knot, cam)
				anchor = S.camera_ray_anchor(cam, ret.x)
				@test anchor.x ≈ ret.x.x atol=1e-12
				@test anchor.y ≈ ret.x.y atol=1e-12
				@test anchor.t ≈ t_obs atol=1e-12
			end
		end

		# FastLight with arbitrary camera direction: simultaneity still holds
		@testset "arbitrary camera direction" begin
			n̂ = normalize(SVector(1.0, 0.0, 1.0))
			cam_arb = S.CameraOrtho(;
				photon_direction=n̂, xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=3.0,
				light=S.FastLight(),
			)
			x0 = S.FourPosition(0.0, 1.0, 0.5, 3.0)
			u = S.FourVelocity(SVector(0.1, 0.05, 0.9))
			wl = S.Geometries.InertialWorldline(x0, u)
			ret = S.retarded_event(wl, cam_arb)
			# FastLight simultaneity: t = t_obs regardless of camera direction
			@test ret.x.t ≈ cam_arb.t atol=1e-12
			# On worldline
			@test ret.x ≈ x0 + u * ret.tau atol=1e-12
			# Roundtrip
			anchor = S.camera_ray_anchor(cam_arb, ret.x)
			@test anchor.t ≈ cam_arb.t atol=1e-12
		end
	end
end

@testitem "RayZ–cone intersection intervals" begin
	import Synchray as S

	ray_at(x, y; ν=1.0) = S.RayZ(; x0=S.FourPosition(0.0, x, y, 0.0), k=ν, nz=16)

	@testset "axis = ẑ (easy hardcode)" begin
		axis = SVector(0.0, 0.0, 1.0)
		φ = 0.2
		s = 1.0..5.0
		
		# On-axis: interval is exactly the truncation.
		@test S._rayz_cone_z_interval(axis, φ, ray_at(0.0, 0.0), s) == s

		# Off-axis but within projected cone: z ≥ x/tanφ.
		x = 0.5 * maximum(s) * tan(φ)
		expected = (x / tan(φ)) .. maximum(s)
		@test S._rayz_cone_z_interval(axis, φ, ray_at(x, 0.0), s) ≈ expected

		# Outside projected cone: miss.
		xmiss = 1.1 * maximum(s) * tan(φ)
		@test S._rayz_cone_z_interval(axis, φ, ray_at(xmiss, 0.0), s) |> isempty
	end

	@testset "tilted axis: z-axis included/excluded" begin
		s = 2.0..6.0

		# If θ < φ, the entire z-axis lies inside the cone (constant opening angle).
		θ_in = 0.05
		φ = 0.2
		axis_in = SVector(sin(θ_in), 0.0, cos(θ_in))
		expected_in = (minimum(s) / cos(θ_in)) .. (maximum(s) / cos(θ_in))
		@test S._rayz_cone_z_interval(axis_in, φ, ray_at(0.0, 0.0), s) ≈ expected_in

		# If θ > φ, the z-axis is outside the cone: miss.
		θ_out = 0.3
		φ = 0.1
		axis_out = SVector(sin(θ_out), 0.0, cos(θ_out))
		@test S._rayz_cone_z_interval(axis_out, φ, ray_at(0.0, 0.0), s) |> isempty
	end

	@testset "axis ⟂ ray direction (axis = x̂)" begin
		axis = SVector(1.0, 0.0, 0.0)
		φ = 0.3
		s = 1.0..3.0

		# For axis=x̂, s(z)=x is constant; require x∈s and y^2+z^2 ≤ tan^2(φ)*x^2.
		x = 2.0
		y = 0.0
		dz = abs(x) * tan(φ)
		expected = (-dz) .. dz
		@test S._rayz_cone_z_interval(axis, φ, ray_at(x, y), s) ≈ expected

		# Same cone but outside truncation in s: miss.
		@test S._rayz_cone_z_interval(axis, φ, ray_at(4.0, 0.0), s) |> isempty

		# av==0 and s0<0: rejected by the +axis half-cone.
		@test S._rayz_cone_z_interval(axis, φ, ray_at(-2.0, 0.0), s) |> isempty
	end

	@testset "edge cases" begin
		# Ray goes exactly along the cone surface (a generatrix): z-axis is on the boundary.
		# Choose axis tilted by θ==φ so that angle(axis, e_z)==φ.
		φ = 0.2
		axis = SVector(sin(φ), 0.0, cos(φ))
		s = 1.0..5.0
		expected = (minimum(s) / cos(φ)) .. (maximum(s) / cos(φ))
		@test S._rayz_cone_z_interval(axis, φ, ray_at(0.0, 0.0), s) ≈ expected

		# Apex-only contact: z-axis outside the cone for θ>φ, but still intersects at z=0.
		# With s including 0, the intersection is a single point.
		θ = 0.3
		φ = 0.1
		axis = SVector(sin(θ), 0.0, cos(θ))
		s = 0.0..1.0
		@test S._rayz_cone_z_interval(axis, φ, ray_at(0.0, 0.0), s) == 0.0 .. 0.0
	end
end


@testitem "Jet point queries helpers" begin
	import Synchray as S

	@testset "event_on_camera_ray matches RayZ convention" begin
		cam = S.CameraZ(; xys=SVector{2}[(0.0, 0.0)], nz=8, ν=1.0, t=2.0)
		r = SVector(0.5, -0.25, 3.0)
		x4 = S.event_on_camera_ray(cam, r)
		@test x4 == S.FourPosition(cam.t + r.z, r.x, r.y, r.z)
	end

	@testset "is_inside_jet (boundary inclusive)" begin
		φj = 0.2
		jet = S.ConicalJet(;
			axis=SVector(0.0, 0.0, 1.0),
			φj,
			s=0.0..10.0,
			ne=S.PowerLaw(-2; val0=1.0, s0=1.0),
			B=S.BFieldSpec_OLD(S.PowerLaw(-1; val0=1.0, s0=1.0), S.ScalarField(), b -> S.FullyTangled(b)),
			speed_profile=(η -> (S.beta, 0.0)),
			model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1.0),
		)

		s = 5.0
		ρmax = s * tan(φj)
		x_inside = S.FourPosition(0.0, 0.5 * ρmax, 0.0, s)
		x_boundary = S.FourPosition(0.0, ρmax, 0.0, s)
		x_outside = S.FourPosition(0.0, 1.01 * ρmax, 0.0, s)
		x_below = S.FourPosition(0.0, 0.0, 0.0, -1.0)

		@test S.is_inside_jet(jet, x_inside)
		@test S.is_inside_jet(jet, x_boundary)
		@test !S.is_inside_jet(jet, x_outside)
		@test !S.is_inside_jet(jet, x_below)

		jetp = S.JetWithPatterns(jet, ())
		@test S.is_inside_jet(jetp, x_inside) == S.is_inside_jet(jet, x_inside)
	end
end


@testitem "Coordinate transform helpers" begin
	import Synchray as S

	@testset "lab <-> jet coords roundtrip" begin
		axis = normalize(SVector(0.2, -0.7, 0.4))
		r = SVector(1.1, -0.3, 2.7)
		rj = S.lab_to_jet_coords(axis, r)
		r2 = S.jet_to_lab_coords(axis, rj)
		@test r2 ≈ r
		@test rj.z ≈ dot(r, axis)
	end

	@testset "jet_rotation_matrix is minimal rotation from lab" begin
		@testset for θ in [0, 0.1, π/2, deg2rad(175)]
			# Axis tilted by θ around lab-y: expect ey == ŷ and (ex, ez) rotated by θ.
			axis = normalize(SVector(sin(θ), 0, cos(θ)))
			R = S.jet_rotation_matrix(axis)
			ex, ey, ez = eachcol(R)
			@test ez ≈ axis
			@test ey ≈ SVector(0, 1, 0)
			@test ex ≈ SVector(cos(θ), 0, -sin(θ))

			# Axis tilted by θ around lab-x: expect ex == x̂ and (ey, ez) rotated by θ.
			axis = normalize(SVector(0, sin(θ), cos(θ)))
			R = S.jet_rotation_matrix(axis)
			ex, ey, ez = eachcol(R)
			@test ez ≈ axis
			@test ex ≈ SVector(1, 0, 0)
			@test ey ≈ SVector(0, cos(θ), -sin(θ))
		end
	end

	@testset "lab <-> camera-ray anchor roundtrip" begin
		x4 = S.FourPosition(7.5, -0.2, 1.1, 3.0)
		x0 = S.camera_ray_anchor(x4)
		@test x0.z == 0
		@test x0.x == x4.x
		@test x0.y == x4.y
		@test x0.t == x4.t - x4.z
	end
end

@testitem "ray_in_jet_coords and camera_band_in_jet_coords" begin
	import Synchray as S
	using RectiGrids

	φj = 0.1
	axis = SVector(sin(0.2), 0.0, cos(0.2))
	jet = S.ConicalJet(;
		axis,
		φj,
		s=(1.0 .. 10.0),
		ne=S.PowerLaw(-2; val0=1.0, s0=1.0),
		B=S.BFieldSpec_OLD(S.PowerLaw(-1; val0=1.0, s0=1.0), S.ScalarField(), b -> S.FullyTangled(b)),
		speed_profile=(η -> (S.beta, 0.9)),
		model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=0.1),
	) |> S.prepare_for_computations

	z_range = 0.0 .. 20.0

	@testset "ray_in_jet_coords" begin
		ray = S.RayZ(; x0=S.FourPosition(0.0, 1.0, 0.0, 0.0), k=1.0, nz=16)
		pts = S.ray_in_jet_coords(ray, jet; z_range)

		@test length(pts) == 2
		@test all(p -> p isa SVector{3}, pts)
		@swiz pts.z
		# Line should have different s values at endpoints (s is z-component in jet frame)
		@test pts[1][3] != pts[2][3]
	end

	@testset "camera_band_in_jet_coords" begin
		cam = S.CameraZ(; xys=grid(SVector, x=range(-1.0, 1.0, 8), y=range(-1.0, 1.0, 8)), nz=16, ν=1.0, t=0.0)
		corners = S.camera_band_in_jet_coords(cam, jet; z_range)

		@test length(corners) == 4
		@test all(c -> c isa SVector{3}, corners)
		# Corners should form a quadrilateral (not all same point)
		@test !allequal(corners)
	end
end

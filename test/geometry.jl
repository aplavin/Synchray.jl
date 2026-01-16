@testitem "Conical geometry" begin
    import Synchray as S

	axis = SVector(0, 0, 1)
	φj = 0.05
	z = 1.0 .. 5.0
	
	geom = Geometries.Conical(axis, φj, z)
	@test S.geometry_axis(geom) == axis
	
	# Test natural_coords
	x4_on_axis = S.FourPosition(0, 0, 0, 2.0)
	coords = S.natural_coords(geom, x4_on_axis)
	@test coords.z ≈ 2.0
	@test coords.ρ ≈ 0.0
	@test coords.η ≈ 0.0
	
	# Test Val(:z) optimization
	@test S.natural_coords(geom, x4_on_axis, Val(:z)) ≈ 2.0
	
	# Off-axis point
	x4_off = S.FourPosition(0, 0.1, 0, 2.0)
	coords_off = S.natural_coords(geom, x4_off)
	@test coords_off.z ≈ 2.0
	@test coords_off.ρ ≈ 0.1
	@test coords_off.η ≈ 0.1 / (2.0 * tan(φj))
	
	# Test is_inside
	@test S.is_inside(geom, x4_on_axis)
	# Point outside s_range
	x4_out = S.FourPosition(0, 0, 0, 10.0)
	@test !S.is_inside(geom, x4_out)
	# Point outside cone but in s_range
	x4_wide = S.FourPosition(0, 1.0, 0, 2.0)
	@test !S.is_inside(geom, x4_wide)
	
	# Test z_interval (on-axis ray)
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=1024)
	@test S.z_interval(geom, ray) == z
	
	# prepare_for_computations should cache trig
	geom_prepared = S.prepare_for_computations(geom)
	@test tan(geom_prepared.φj) ≈ tan(φj)
	@test cos(geom_prepared.φj) ≈ cos(φj)
	
	# Test Accessors.set for geometry_axis
	new_axis = normalize(SVector(1, 0, 1))
	geom_new = @set S.geometry_axis(geom) = new_axis
	@test S.geometry_axis(geom_new) == new_axis
	@test S.geometry_axis(geom) == axis  # original unchanged
end

@testitem "Coordinate transformations" begin
    import Synchray as S

	# Basic roundtrip test
	axis = normalize(SVector(1, 0, 1))
	geom = S.Geometries.Conical(; axis, φj=0.1, z=1.0 .. 5.0)
	
	# Test rotation matrix
	R = S.rotation_local_to_lab(geom)
	@test size(R) == (3, 3)
	@test det(R) ≈ 1  # proper rotation
	
	# Test roundtrip
	r_lab = SVector(1.0, 2.0, 3.0)
	r_local = S.rotate_lab_to_local(geom, r_lab)
	r_back = S.rotate_local_to_lab(geom, r_local)
	@test r_back ≈ r_lab
	
	# For axis-aligned geometry, local z should align with axis
	geom_z = S.Geometries.Conical(; axis=SVector(0, 0, 1), φj=0.1, z=1.0 .. 5.0)
	R_z = S.rotation_local_to_lab(geom_z)
	@test R_z[:, 3] ≈ SVector(0, 0, 1)  # ez = axis
	
	# Arbitrary axis roundtrip
	axis_arb = normalize(SVector(0.2, -0.7, 0.4))
	geom_arb = Geometries.Conical(axis_arb, 0.1, 1.0 .. 5.0)
	r_arb = SVector(1.1, -0.3, 2.7)
	r_local_arb = S.rotate_lab_to_local(geom_arb, r_arb)
	r_back_arb = S.rotate_local_to_lab(geom_arb, r_local_arb)
	@test r_back_arb ≈ r_arb
	@test r_local_arb[3] ≈ dot(r_arb, axis_arb)  # local z = axial component
	
	# Rotation matrix is minimal rotation from lab
	# Axis tilted by θ around lab-y: expect ey == ŷ and (ex, ez) rotated by θ
	@testset for θ in [0, 0.1, π/2, deg2rad(175)]
		axis_y = normalize(SVector(sin(θ), 0, cos(θ)))
		geom_y = Geometries.Conical(axis_y, 0.1, 1.0 .. 5.0)
		R_y = S.rotation_local_to_lab(geom_y)
		ex, ey, ez = eachcol(R_y)
		@test ez ≈ axis_y
		@test ey ≈ SVector(0, 1, 0)
		@test ex ≈ SVector(cos(θ), 0, -sin(θ))
		
		# Axis tilted by θ around lab-x: expect ex == x̂ and (ey, ez) rotated by θ
		axis_x = normalize(SVector(0, sin(θ), cos(θ)))
		geom_x = Geometries.Conical(axis_x, 0.1, 1.0 .. 5.0)
		R_x = S.rotation_local_to_lab(geom_x)
		ex, ey, ez = eachcol(R_x)
		@test ez ≈ axis_x
		@test ex ≈ SVector(1, 0, 0)
		@test ey ≈ SVector(0, cos(θ), -sin(θ))
	end
end

@testitem "Visualization helpers" begin
    import Synchray as S
    using RectiGrids
	
	φj = 0.1
	axis = SVector(sin(0.2), 0.0, cos(0.2))
	geom = Geometries.Conical(axis, φj, 1.0 .. 10.0)
	z_range = 0.0 .. 20.0
	
	# Test ray_in_local_coords
	ray = S.RayZ(; x0=S.FourPosition(0.0, 1.0, 0.0, 0.0), k=1.0, nz=16)
	pts = S.ray_in_local_coords(ray, geom; z_range)
	
	@test length(pts) == 2
	@test all(p -> p isa SVector{3}, pts)
	# Line should have different s values at endpoints (s is z-component in local frame)
	@test pts[1][3] != pts[2][3]
	
	# Test camera_fov_in_local_coords
	cam = S.CameraZ(; xys=grid(SVector, x=range(-1.0, 1.0, 8), y=range(-1.0, 1.0, 8)), nz=16, ν=1.0, t=0.0)
	corners = S.camera_fov_in_local_coords(cam, geom; z_range)
	
	@test length(corners) == 4
	@test all(c -> c isa SVector{3}, corners)
	# Corners should form a quadrilateral (not all same point)
	@test !allequal(corners)
end

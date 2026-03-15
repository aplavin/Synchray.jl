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
	s_range = 0.0 .. 20.0

	# Test ray_in_local_coords
	ray = S.RayZ(; x0=S.FourPosition(0.0, 1.0, 0.0, 0.0), k=1.0, nz=16)
	pts = S.ray_in_local_coords(ray, geom; s_range)

	@test length(pts) == 2
	@test all(p -> p isa SVector{3}, pts)
	# Line should have different s values at endpoints (s is z-component in local frame)
	@test pts[1][3] != pts[2][3]

	# Test camera_fov_in_local_coords
	cam = S.CameraZ(; xys=grid(SVector, x=range(-1.0, 1.0, 8), y=range(-1.0, 1.0, 8)), nz=16, ν=1.0, t=0.0)
	corners = S.camera_fov_in_local_coords(cam, geom; s_range)
	
	@test length(corners) == 4
	@test all(c -> c isa SVector{3}, corners)
	# Corners should form a quadrilateral (not all same point)
	@test !allequal(corners)
end

@testitem "Arbitrary-angle Ray and Camera basics" begin
	import Synchray as S

	x0 = S.FourPosition(0.0, 0.0, 0.0, 0.0)
	k = S.photon_k(2.0, SVector(0.0, 0.0, 1.0))

	# GPU compatibility: Ray must be isbits
	ray = S.Ray(x0, k, SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), 64)
	@test isbits(ray)

	# RayZ produces a Ray
	rayz = S.RayZ(x0, k, 64)
	@test rayz isa S.Ray
	@test isbits(rayz)

	# direction3 extracts the unit spatial direction
	@test S.direction3(ray) ≈ SVector(0.0, 0.0, 1.0)

	# Arbitrary direction
	n̂ = normalize(SVector(1.0, 0.0, 1.0))
	k_arb = S.photon_k(2.0, n̂)
	e1 = normalize(cross(SVector(0.0, 1.0, 0.0), n̂))
	e2 = cross(n̂, e1)
	ray_arb = S.Ray(x0, k_arb, e1, e2, 64)
	@test isbits(ray_arb)
	@test S.direction3(ray_arb) ≈ n̂
	@test S.frequency(ray_arb) ≈ 2.0

	# Camera with arbitrary look_direction: orthonormal screen basis
	cam = S.CameraOrtho(;
		look_direction=SVector(1.0, 0.0, 1.0),
		xys=SVector{2}[(0.0, 0.0)],
		nz=64, ν=2.0, t=0.0,
	)
	@test cam.n ≈ n̂
	@test norm(cam.e1) ≈ 1
	@test norm(cam.e2) ≈ 1
	@test dot(cam.e1, cam.n) ≈ 0 atol=√eps(1.0)
	@test dot(cam.e2, cam.n) ≈ 0 atol=√eps(1.0)
	@test dot(cam.e1, cam.e2) ≈ 0 atol=√eps(1.0)
	@test cross(cam.n, cam.e1) ≈ cam.e2
end

@testitem "Ellipsoid geometry" begin
	import Synchray as S

	u0 = S.FourVelocity(1.0, 0.0, 0.0, 0.0)  # stationary
	wl = S.Geometries.InertialWorldline(S.FourPosition(0.0, 0.0, 0.0, 0.0), u0)
	sizes = SVector(1.0, 1.0, 1.0)
	geom = S.Geometries.Ellipsoid(wl, sizes)

	# four_velocity returns the worldline velocity
	@test S.four_velocity(geom, S.FourPosition(0,0,0,0)) === u0

	# prepare_for_computations is identity
	@test S.prepare_for_computations(geom) === geom

	# z_interval: ray through center of stationary unit sphere
	ray_center = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=1.0, nz=64)
	zi = S.z_interval(geom, ray_center)
	@test leftendpoint(zi) ≈ -1.0
	@test rightendpoint(zi) ≈ 1.0

	# z_interval: non-spherical sizes — z semi-axis = 2
	geom_elong = S.Geometries.Ellipsoid(wl, SVector(1.0, 1.0, 2.0))
	zi_elong = S.z_interval(geom_elong, ray_center)
	@test leftendpoint(zi_elong) ≈ -2.0
	@test rightendpoint(zi_elong) ≈ 2.0

	# z_interval: ray that misses
	ray_miss = S.RayZ(; x0=S.FourPosition(0.0, 5.0, 0.0, 0.0), k=1.0, nz=64)
	zi_miss = S.z_interval(geom, ray_miss)
	@test isempty(zi_miss)

	# z_interval: moving ellipsoid — shifted center
	wl_shifted = S.Geometries.InertialWorldline(S.FourPosition(0.0, 0.0, 0.0, 3.0), u0)
	geom_shifted = S.Geometries.Ellipsoid(wl_shifted, sizes)
	ray_at_shifted = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=1.0, nz=64)
	zi_shifted = S.z_interval(geom_shifted, ray_at_shifted)
	@test leftendpoint(zi_shifted) ≈ 2.0
	@test rightendpoint(zi_shifted) ≈ 4.0

	# z_interval: moving ellipsoid with nonzero velocity
	β = SVector(0.0, 0.0, 0.5)
	u_moving = S.FourVelocity(β)
	wl_moving = S.Geometries.InertialWorldline(S.FourPosition(0.0, 0.0, 0.0, 0.0), u_moving)
	geom_moving = S.Geometries.Ellipsoid(wl_moving, sizes)
	zi_moving = S.z_interval(geom_moving, ray_center)
	@test !isempty(zi_moving)
	@test S.four_velocity(geom_moving, S.FourPosition(0,0,0,0)) === u_moving
end

@testitem "SlowLight/FastLight basics" begin
	import Synchray as S
	using Accessors
	using RectiGrids

	x0 = S.FourPosition(0.0, 0.0, 0.0, 0.0)

	# isbits (GPU-safe)
	@test isbits(S.SlowLight())
	@test isbits(S.FastLight())
	ray_slow = S.RayZ(; x0, k=2.0, nz=64)
	ray_fast = S.RayZ(; x0, k=2.0, nz=64, light=S.FastLight())
	@test isbits(ray_slow)
	@test isbits(ray_fast)

	# Constructors default to SlowLight
	@test ray_slow.light === S.SlowLight()
	@test S.CameraZ(; xys=grid(SVector, x=[0.0], y=[0.0]), nz=1, ν=1.0, t=0.0).light === S.SlowLight()
	@test S.CameraOrtho(; look_direction=SVector(1.0, 0.0, 1.0), xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=0.0).light === S.SlowLight()
	k = S.photon_k(2.0, SVector(0.0, 0.0, 1.0))
	@test S.Ray(x0, k, SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), 64).light === S.SlowLight()

	# direction4(ray): Z-ray
	d_slow = S.direction4(ray_slow)
	d_fast = S.direction4(ray_fast)
	@test d_slow.t == 1.0
	@test d_fast.t == 0.0
	@test SVector(d_slow.x, d_slow.y, d_slow.z) ≈ S.direction3(ray_slow)
	@test SVector(d_fast.x, d_fast.y, d_fast.z) ≈ S.direction3(ray_fast)

	# direction4(ray): arbitrary directions
	@testset for n̂ in [normalize(SVector(1.0, 0.0, 0.0)),
	                    normalize(SVector(1.0, 1.0, 1.0)),
	                    normalize(SVector(-0.3, 0.7, 0.5))]
		k_arb = S.photon_k(2.0, n̂)
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		e1 = normalize(cross(up, n̂))
		e2 = cross(n̂, e1)
		r_slow = S.Ray(x0, k_arb, e1, e2, 64, S.SlowLight())
		r_fast = S.Ray(x0, k_arb, e1, e2, 64, S.FastLight())
		@test S.direction4(r_slow).t == 1.0
		@test S.direction4(r_fast).t == 0.0
		@test SVector(S.direction4(r_slow).x, S.direction4(r_slow).y, S.direction4(r_slow).z) ≈ n̂
		@test SVector(S.direction4(r_fast).x, S.direction4(r_fast).y, S.direction4(r_fast).z) ≈ n̂
	end

	# event_on_camera_ray: exact formulas for both modes
	t_obs = 2.5
	r = SVector(3.0, 0.5, -0.25)
	@testset for n̂ in [SVector(0.0, 0.0, 1.0),
	                    SVector(1.0, 0.0, 0.0),
	                    normalize(SVector(1.0, 1.0, 1.0)),
	                    normalize(SVector(-0.3, 0.7, 0.5))]
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		cam_slow = S.CameraOrtho(; look_direction=n̂, up, xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=t_obs)
		cam_fast = S.CameraOrtho(; look_direction=n̂, up, xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=t_obs, light=S.FastLight())

		s = dot(r - cam_slow.origin, cam_slow.n)

		x4_slow = S.event_on_camera_ray(cam_slow, r)
		x4_fast = S.event_on_camera_ray(cam_fast, r)

		# SlowLight: t = t_obs + depth
		@test x4_slow.t ≈ t_obs + s
		# FastLight: t = t_obs
		@test x4_fast.t ≈ t_obs
		# Both: same spatial coordinates
		@test SVector(x4_slow.x, x4_slow.y, x4_slow.z) ≈ r
		@test SVector(x4_fast.x, x4_fast.y, x4_fast.z) ≈ r
	end

	# camera_ray_anchor: exact formulas for both modes
	x4 = S.FourPosition(7.5, -0.2, 1.1, 3.0)
	@testset for n̂ in [SVector(0.0, 0.0, 1.0),
	                    SVector(1.0, 0.0, 0.0),
	                    normalize(SVector(1.0, 1.0, 1.0))]
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		cam_slow = S.CameraOrtho(; look_direction=n̂, up, xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=t_obs)
		cam_fast = S.CameraOrtho(; look_direction=n̂, up, xys=SVector{2}[(0.0, 0.0)], nz=1, ν=1.0, t=t_obs, light=S.FastLight())

		r_lab = SVector(x4.x, x4.y, x4.z)
		s = dot(r_lab - cam_slow.origin, cam_slow.n)
		r_screen = r_lab - s * cam_slow.n

		anchor_slow = S.camera_ray_anchor(cam_slow, x4)
		anchor_fast = S.camera_ray_anchor(cam_fast, x4)

		# SlowLight: t_cam = t - depth
		@test anchor_slow.t ≈ x4.t - s
		# FastLight: t_cam = t
		@test anchor_fast.t ≈ x4.t
		# Both: same screen projection
		@test SVector(anchor_slow.x, anchor_slow.y, anchor_slow.z) ≈ r_screen
		@test SVector(anchor_fast.x, anchor_fast.y, anchor_fast.z) ≈ r_screen

		# Roundtrip: camera_ray_anchor(event_on_camera_ray(r)).t ≈ t_obs
		@test S.camera_ray_anchor(cam_slow, S.event_on_camera_ray(cam_slow, r_lab)).t ≈ t_obs
		@test S.camera_ray_anchor(cam_fast, S.event_on_camera_ray(cam_fast, r_lab)).t ≈ t_obs
	end
end

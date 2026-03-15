@testitem "z_interval for conical geometry" begin
	import Synchray as S

	ray_at(x, y; ν=1.0) = S.RayZ(; x0=S.FourPosition(0.0, x, y, 0.0), k=ν, nz=16)

	@testset "axis = ẑ (easy hardcode)" begin
		axis = SVector(0.0, 0.0, 1.0)
		φ = 0.2
		z = 1.0..5.0
		geom = S.Geometries.Conical(; axis, φj=φ, z)
		
		# On-axis: interval is exactly the truncation.
		@test S.z_interval(geom, ray_at(0.0, 0.0)) == z

		# Off-axis but within projected cone.
		x = 0.5 * maximum(z) * tan(φ)
		expected = (x / tan(φ)) .. maximum(z)
		@test S.z_interval(geom, ray_at(x, 0.0)) ≈ expected

		# Outside projected cone: miss.
		xmiss = 1.1 * maximum(z) * tan(φ)
		@test S.z_interval(geom, ray_at(xmiss, 0.0)) |> isempty
	end

	@testset "tilted axis: z-axis included/excluded" begin
		z = 2.0..6.0

		# If θ < φ, the entire z-axis lies inside the cone (constant opening angle).
		θ_in = 0.05
		φ = 0.2
		axis_in = SVector(sin(θ_in), 0.0, cos(θ_in))
		geom_in = S.Geometries.Conical(; axis=axis_in, φj=φ, z)
		expected_in = (minimum(z) / cos(θ_in)) .. (maximum(z) / cos(θ_in))
		@test S.z_interval(geom_in, ray_at(0.0, 0.0)) ≈ expected_in

		# If θ > φ, the z-axis is outside the cone: miss.
		θ_out = 0.3
		φ = 0.1
		axis_out = SVector(sin(θ_out), 0.0, cos(θ_out))
		geom_out = S.Geometries.Conical(; axis=axis_out, φj=φ, z)
		@test S.z_interval(geom_out, ray_at(0.0, 0.0)) |> isempty
	end

	@testset "axis ⟂ ray direction (axis = x̂)" begin
		axis = SVector(1.0, 0.0, 0.0)
		φ = 0.3
		z = 1.0..3.0
		geom = S.Geometries.Conical(; axis, φj=φ, z)

		# For axis=x̂, z-coordinate is x (constant); require x∈z and y^2+z^2 ≤ tan^2(φ)*x^2.
		x = 2.0
		y = 0.0
		dz = abs(x) * tan(φ)
		expected = (-dz) .. dz
		@test S.z_interval(geom, ray_at(x, y)) ≈ expected

		# Same cone but outside truncation in z: miss.
		@test S.z_interval(geom, ray_at(4.0, 0.0)) |> isempty

		# Point outside the +axis half-cone.
		@test S.z_interval(geom, ray_at(-2.0, 0.0)) |> isempty
	end

	@testset "edge cases" begin
		# Ray goes exactly along the cone surface (a generatrix): z-axis on the boundary.
		# Choose axis tilted by θ==φ so that angle(axis, e_z)==φ.
		φ = 0.2
		axis = SVector(sin(φ), 0.0, cos(φ))
		z = 1.0..5.0
		geom = S.Geometries.Conical(; axis, φj=φ, z)
		expected = (minimum(z) / cos(φ)) .. (maximum(z) / cos(φ))
		@test S.z_interval(geom, ray_at(0.0, 0.0)) ≈ expected

		# Apex-only contact: z-axis outside the cone for θ>φ, but still intersects at z=0.
		# With z including 0, the intersection is a single point.
		θ = 0.3
		φ = 0.1
		axis = SVector(sin(θ), 0.0, cos(θ))
		z = 0.0..1.0
		geom = S.Geometries.Conical(; axis, φj=φ, z)
		@test S.z_interval(geom, ray_at(0.0, 0.0)) == 0.0 .. 0.0
	end
end


@testitem "z_interval for conical geometry: non-Z rays" begin
	import Synchray as S

	# Helper: create a Ray at position r0 with direction n̂
	function arb_ray(r0, n̂; ν=1.0, nz=16)
		n̂ = normalize(n̂)
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		e1 = normalize(cross(up, n̂))
		e2 = cross(n̂, e1)
		S.Ray(S.FourPosition(0.0, r0...), S.photon_k(ν, n̂), e1, e2, nz)
	end

	@testset "Tilted ray through apex of Z-cone" begin
		φj = 0.2
		z = 1.0 .. 5.0
		geom = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj, z)

		# Ray from origin at angle θ < φj: inside cone, interval = [z1/cos(θ), z2/cos(θ)]
		θ_in = 0.05
		n̂_in = SVector(sin(θ_in), 0.0, cos(θ_in))
		@test S.z_interval(geom, arb_ray(SA[0,0,0], n̂_in)) ≈ (1.0 / cos(θ_in)) .. (5.0 / cos(θ_in))

		# Ray from origin at angle θ > φj: misses
		θ_out = 0.3
		n̂_out = SVector(sin(θ_out), 0.0, cos(θ_out))
		@test S.z_interval(geom, arb_ray(SA[0,0,0], n̂_out)) |> isempty

		# Ray at exactly θ = φj: on the cone surface (boundary)
		n̂_edge = SVector(sin(φj), 0.0, cos(φj))
		@test S.z_interval(geom, arb_ray(SA[0,0,0], n̂_edge)) ≈ (1.0 / cos(φj)) .. (5.0 / cos(φj))
	end

	@testset "Horizontal ray through Z-cone" begin
		φj = 0.3
		z = 1.0 .. 5.0
		geom = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj, z)

		# Horizontal ray (along X) at height z_mid: crosses cone where |x| ≤ z_mid * tan(φj)
		z_mid = 3.0
		r_cone = z_mid * tan(φj)
		@test S.z_interval(geom, arb_ray(SA[0,0,z_mid], SA[1,0,0])) ≈ (-r_cone) .. r_cone

		# Horizontal ray above truncation: misses
		@test S.z_interval(geom, arb_ray(SA[0,0,6.0], SA[1,0,0])) |> isempty

		# Horizontal ray below truncation: misses
		@test S.z_interval(geom, arb_ray(SA[0,0,0.5], SA[1,0,0])) |> isempty
	end

	@testset "45° ray through Z-cone" begin
		z = 1.0 .. 5.0
		n̂_45 = normalize(SVector(1.0, 0.0, 1.0))

		# For φj < π/4: the 45° ray is outside the cone (tan(45°) = 1 > tan(φj))
		geom_narrow = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj=0.4, z)
		@test S.z_interval(geom_narrow, arb_ray(SA[0,0,0], n̂_45)) |> isempty

		# For φj > π/4: the 45° ray is inside the cone
		geom_wide = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj=1.0, z)
		@test S.z_interval(geom_wide, arb_ray(SA[0,0,0], n̂_45)) ≈ (1.0 * √2) .. (5.0 * √2)
	end

	@testset "Tilted cone with tilted ray: consistency with Z-ray equivalent" begin
		# A tilted cone seen by a Z-ray should give the same result
		# as a Z-cone seen by a tilted ray (after coordinate rotation).
		φj = 0.2
		z = 1.0 .. 5.0
		θ = 0.05  # tilt angle

		# Case 1: tilted cone, Z-ray from origin
		axis_tilted = SVector(sin(θ), 0.0, cos(θ))
		geom_tilted = S.Geometries.Conical(; axis=axis_tilted, φj, z)
		z_ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=1.0, nz=16)
		seg1 = S.z_interval(geom_tilted, z_ray)

		# Case 2: Z-cone, tilted ray from origin (same tilt angle)
		geom_z = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj, z)
		n̂_tilted = SVector(sin(θ), 0.0, cos(θ))
		seg2 = S.z_interval(geom_z, arb_ray(SA[0,0,0], n̂_tilted))

		# Both should give the same interval width (same geometry, just rotated)
		@test S.width(seg1) ≈ S.width(seg2)
	end

	@testset "event_on_camera_ray with arbitrary Camera" begin
		cam = S.Camera(;
			look_direction=SVector(1.0, 0.0, 0.0),
			origin=SVector(0.0, 0.0, 0.0),
			xys=SVector{2}[(0.0, 0.0)], nz=8, ν=1.0, t=2.0,
		)
		r = SVector(3.0, 0.5, -0.25)
		x4 = S.event_on_camera_ray(cam, r)
		# s = dot(r - origin, n) = dot(r, x̂) = r.x = 3
		@test x4.t ≈ cam.t + 3.0
		@test @swiz(x4.xyz) ≈ r

		# camera_ray_anchor roundtrip
		x4_test = S.FourPosition(7.5, 3.0, -0.2, 1.1)
		x0 = S.camera_ray_anchor(cam, x4_test)
		# s = dot(r - origin, n) = r.x = 3.0
		# r_screen = origin + (r - origin) - s * n = r - 3*x̂ = (0, -0.2, 1.1)
		@test x0.x ≈ 0 atol=1e-14  # screen-plane x (along n̂) should be zero
		@test x0.y ≈ x4_test.y
		@test x0.z ≈ x4_test.z
		@test x0.t ≈ x4_test.t - 3.0
	end
end


@testitem "event_on_camera_ray matches RayZ convention" begin
	import Synchray as S
	
	cam = S.CameraZ(; xys=SVector{2}[(0.0, 0.0)], nz=8, ν=1.0, t=2.0)
	r = SVector(0.5, -0.25, 3.0)
	x4 = S.event_on_camera_ray(cam, r)
	@test x4 == S.FourPosition(cam.t + r.z, r.x, r.y, r.z)
end


@testitem "lab <-> camera-ray anchor roundtrip" begin
	import Synchray as S
	using RectiGrids

	cam = S.CameraZ(; xys=grid(SVector, x=[0.0], y=[0.0]), nz=1, ν=1.0, t=0.0)
	x4 = S.FourPosition(7.5, -0.2, 1.1, 3.0)
	x0 = S.camera_ray_anchor(cam, x4)
	@test x0.z ≈ 0 atol=1e-14
	@test x0.x == x4.x
	@test x0.y == x4.y
	@test x0.t == x4.t - x4.z
end

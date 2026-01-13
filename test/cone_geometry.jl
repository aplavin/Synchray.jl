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


@testitem "event_on_camera_ray matches RayZ convention" begin
	import Synchray as S
	
	cam = S.CameraZ(; xys=SVector{2}[(0.0, 0.0)], nz=8, ν=1.0, t=2.0)
	r = SVector(0.5, -0.25, 3.0)
	x4 = S.event_on_camera_ray(cam, r)
	@test x4 == S.FourPosition(cam.t + r.z, r.x, r.y, r.z)
end


@testitem "lab <-> camera-ray anchor roundtrip" begin
	import Synchray as S

	x4 = S.FourPosition(7.5, -0.2, 1.1, 3.0)
	x0 = S.camera_ray_anchor(x4)
	@test x0.z == 0
	@test x0.x == x4.x
	@test x0.y == x4.y
	@test x0.t == x4.t - x4.z
end

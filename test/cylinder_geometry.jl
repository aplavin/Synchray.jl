@testitem "Cylindrical geometry" begin
	import Synchray as S
	using Accessors

	axis = SVector(0, 0, 1)
	radius = 0.5
	z = 1.0 .. 5.0

	geom = S.Geometries.Cylindrical(; axis, radius, z)
	@test S.geometry_axis(geom) == axis

	# natural_coords on-axis
	x4_on_axis = S.FourPosition(0, 0, 0, 2.0)
	coords = S.natural_coords(geom, x4_on_axis)
	@test coords.z ≈ 2.0
	@test coords.ρ ≈ 0.0
	@test coords.η ≈ 0.0

	# Val(:z) optimization
	@test S.natural_coords(geom, x4_on_axis, Val(:z)) ≈ 2.0

	# Off-axis point: η = ρ / radius (constant, unlike cone's z-dependent η)
	x4_off = S.FourPosition(0, 0.1, 0, 2.0)
	coords_off = S.natural_coords(geom, x4_off)
	@test coords_off.z ≈ 2.0
	@test coords_off.ρ ≈ 0.1
	@test coords_off.η ≈ 0.1 / radius

	# Same off-axis point at different z: η should be the same (cylinder, not cone)
	x4_off2 = S.FourPosition(0, 0.1, 0, 4.0)
	@test S.natural_coords(geom, x4_off2).η ≈ coords_off.η

	# is_inside
	@test S.is_inside(geom, x4_on_axis)
	# Outside z range
	@test !S.is_inside(geom, S.FourPosition(0, 0, 0, 10.0))
	# Outside radius but in z range
	@test !S.is_inside(geom, S.FourPosition(0, 1.0, 0, 2.0))
	# At radius boundary
	@test S.is_inside(geom, S.FourPosition(0, radius, 0, 2.0))

	# z_interval on-axis ray
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=1024)
	@test S.z_interval(geom, ray) == z

	# prepare_for_computations is identity for Cylindrical
	@test S.prepare_for_computations(geom) === geom

	# Accessors.set for geometry_axis
	new_axis = normalize(SVector(1, 0, 1))
	geom_new = @set S.geometry_axis(geom) = new_axis
	@test S.geometry_axis(geom_new) ≈ new_axis
	@test S.geometry_axis(geom) == axis  # original unchanged
	@test geom_new.radius == radius  # preserved
	@test geom_new.z == z  # preserved

	# rotation_local_to_lab
	R = S.rotation_local_to_lab(geom)
	@test size(R) == (3, 3)
	@test det(R) ≈ 1
end


@testitem "z_interval for cylindrical geometry" begin
	import Synchray as S

	ray_at(x, y; ν=1.0) = S.RayZ(; x0=S.FourPosition(0.0, x, y, 0.0), k=ν, nz=16)

	@testset "axis = ẑ (easy hardcode)" begin
		axis = SVector(0.0, 0.0, 1.0)
		R = 1.0
		z = 1.0..5.0
		geom = S.Geometries.Cylindrical(; axis, radius=R, z)

		# On-axis: interval is exactly the truncation
		@test S.z_interval(geom, ray_at(0.0, 0.0)) == z

		# Off-axis but within cylinder: still full z range (constant radius)
		@test S.z_interval(geom, ray_at(0.5, 0.0)) == z

		# At cylinder boundary: parallel ray touching surface
		@test S.z_interval(geom, ray_at(R, 0.0)) ≈ z

		# Outside cylinder: miss
		@test S.z_interval(geom, ray_at(1.1 * R, 0.0)) |> isempty

		# Off-axis in both x and y, within cylinder
		@test S.z_interval(geom, ray_at(0.5, 0.5)) == z

		# Off-axis in both x and y, outside cylinder
		@test S.z_interval(geom, ray_at(0.8, 0.8)) |> isempty
	end

	@testset "tilted axis" begin
		z = 2.0..6.0
		R = 1.0

		# Small tilt: z-axis is inside the cylinder (perpendicular distance = 0)
		θ_in = 0.05
		axis_in = SVector(sin(θ_in), 0.0, cos(θ_in))
		geom_in = S.Geometries.Cylindrical(; axis=axis_in, radius=R, z)
		expected_in = (minimum(z) / cos(θ_in)) .. (maximum(z) / cos(θ_in))
		@test S.z_interval(geom_in, ray_at(0.0, 0.0)) ≈ expected_in

		# Axis perpendicular to z: axis=x̂, z-ray at origin has axial coord x=0,
		# which is outside z=2..6 → miss
		axis_perp = SVector(1.0, 0.0, 0.0)
		geom_perp = S.Geometries.Cylindrical(; axis=axis_perp, radius=R, z)
		@test S.z_interval(geom_perp, ray_at(0.0, 0.0)) |> isempty
	end

	@testset "axis ⟂ ray direction (axis = x̂)" begin
		axis = SVector(1.0, 0.0, 0.0)
		R = 0.5
		z = 1.0..3.0
		geom = S.Geometries.Cylindrical(; axis, radius=R, z)

		# Z-ray at (x=2, y=0): axial coord = x = 2 ∈ [1,3].
		# r0_perp = (0, 0, 0), n̂_perp = (0, 0, 1).
		# A=1, B=0, C=-R². s ∈ [-R, R].
		@test S.z_interval(geom, ray_at(2.0, 0.0)) ≈ (-R) .. R

		# Outside truncation: miss
		@test S.z_interval(geom, ray_at(4.0, 0.0)) |> isempty

		# Negative axis direction: outside positive z range
		@test S.z_interval(geom, ray_at(-2.0, 0.0)) |> isempty

		# Z-ray at (x=2, y=0.3): r0_perp = (0, 0.3, 0).
		# C = 0.3² - 0.5² = -0.16. D = 4*0.16 = 0.64.
		# s ∈ [-√0.16, √0.16] = [-0.4, 0.4]
		@test S.z_interval(geom, ray_at(2.0, 0.3)) ≈ (-0.4) .. 0.4
	end

	@testset "key difference from cone: z-independence" begin
		axis = SVector(0.0, 0.0, 1.0)
		R = 1.0
		z = 1.0..10.0
		geom = S.Geometries.Cylindrical(; axis, radius=R, z)

		# For a cylinder, any Z-ray within radius hits the full z range,
		# regardless of where we are. This is unlike a cone where the
		# entry point depends on the cone opening angle.
		@test S.z_interval(geom, ray_at(0.3, 0.0)) == z
		@test S.z_interval(geom, ray_at(0.9, 0.0)) == z
		@test S.z_interval(geom, ray_at(0.0, 0.9)) == z
	end
end


@testitem "z_interval for cylindrical geometry: non-Z rays" begin
	import Synchray as S

	function arb_ray(r0, n̂; ν=1.0, nz=16)
		n̂ = normalize(n̂)
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		e1 = normalize(cross(up, n̂))
		e2 = cross(n̂, e1)
		S.Ray(S.FourPosition(0.0, r0...), S.photon_k(ν, n̂), e1, e2, nz)
	end

	@testset "Horizontal ray through Z-cylinder" begin
		R = 1.0
		z = 1.0 .. 5.0
		geom = S.Geometries.Cylindrical(; axis=SVector(0.0, 0.0, 1.0), radius=R, z)

		# Horizontal ray (along X) at height z_mid: chord = 2R
		z_mid = 3.0
		@test S.z_interval(geom, arb_ray(SA[0,0,z_mid], SA[1,0,0])) ≈ (-R) .. R

		# Above truncation: miss
		@test S.z_interval(geom, arb_ray(SA[0,0,6.0], SA[1,0,0])) |> isempty

		# Below truncation: miss
		@test S.z_interval(geom, arb_ray(SA[0,0,0.5], SA[1,0,0])) |> isempty

		# Offset in y: chord = 2√(R²-y²)
		y_off = 0.5
		chord = √(R^2 - y_off^2)
		@test S.z_interval(geom, arb_ray(SA[0,y_off,z_mid], SA[1,0,0])) ≈ (-chord) .. chord

		# Tangent horizontal ray at y=R: D=0 → empty (consistent with conical behavior)
		@test S.z_interval(geom, arb_ray(SA[0,R,z_mid], SA[1,0,0])) |> isempty

		# Horizontal ray outside cylinder: miss
		@test S.z_interval(geom, arb_ray(SA[0,1.5,z_mid], SA[1,0,0])) |> isempty
	end

	@testset "Tilted ray through Z-cylinder" begin
		R = 1.0
		z = 1.0 .. 5.0
		geom = S.Geometries.Cylindrical(; axis=SVector(0.0, 0.0, 1.0), radius=R, z)

		# Ray from origin along axis: full z range
		@test S.z_interval(geom, arb_ray(SA[0,0,0], SA[0,0,1])) ≈ z

		# Ray from origin at angle θ: perpendicular distance ρ(s) = s·sin(θ),
		# so cylinder constraint: s ≤ R/sin(θ).
		# Axial truncation: s·cos(θ) ∈ [1,5], so s ∈ [1/cos(θ), 5/cos(θ)].
		θ = 0.1
		n̂ = SVector(sin(θ), 0.0, cos(θ))
		s_trunc_lo = 1.0 / cos(θ)
		s_trunc_hi = 5.0 / cos(θ)
		s_cyl = R / sin(θ)
		expected = s_trunc_lo .. min(s_trunc_hi, s_cyl)
		@test S.z_interval(geom, arb_ray(SA[0,0,0], n̂)) ≈ expected
	end

	@testset "Consistency: tilted cylinder + Z-ray vs Z-cylinder + tilted ray" begin
		R = 1.0
		z = 1.0 .. 5.0
		θ = 0.05

		# Case 1: tilted cylinder, Z-ray from origin
		axis_tilted = SVector(sin(θ), 0.0, cos(θ))
		geom_tilted = S.Geometries.Cylindrical(; axis=axis_tilted, radius=R, z)
		z_ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=1.0, nz=16)
		seg1 = S.z_interval(geom_tilted, z_ray)

		# Case 2: Z-cylinder, tilted ray from origin
		geom_z = S.Geometries.Cylindrical(; axis=SVector(0.0, 0.0, 1.0), radius=R, z)
		n̂_tilted = SVector(sin(θ), 0.0, cos(θ))
		seg2 = S.z_interval(geom_z, arb_ray(SA[0,0,0], n̂_tilted))

		# Same geometry, just rotated → same interval width
		@test S.width(seg1) ≈ S.width(seg2)
	end

	@testset "Cylinder vs cone comparison at matching cross-section" begin
		# A horizontal ray through both at z=z_mid should give the same chord
		# when the cone's radius at z_mid equals the cylinder's radius.
		R = 1.0
		z = 1.0 .. 5.0
		z_mid = 3.0

		geom_cyl = S.Geometries.Cylindrical(; axis=SVector(0.0, 0.0, 1.0), radius=R, z)
		zi_cyl = S.z_interval(geom_cyl, arb_ray(SA[0,0,z_mid], SA[1,0,0]))
		@test zi_cyl ≈ (-R) .. R

		# Cone with tan(φj) = R/z_mid → same radius at z_mid
		φj = atan(R / z_mid)
		geom_cone = S.Geometries.Conical(; axis=SVector(0.0, 0.0, 1.0), φj, z)
		zi_cone = S.z_interval(geom_cone, arb_ray(SA[0,0,z_mid], SA[1,0,0]))
		@test zi_cone ≈ (-R) .. R
		@test zi_cyl ≈ zi_cone
	end
end


@testitem "Cylindrical coordinate transformations" begin
	import Synchray as S

	# Roundtrip test with tilted axis
	axis = normalize(SVector(1, 0, 1))
	geom = S.Geometries.Cylindrical(; axis, radius=0.5, z=1.0 .. 5.0)

	R = S.rotation_local_to_lab(geom)
	@test size(R) == (3, 3)
	@test det(R) ≈ 1  # proper rotation

	r_lab = SVector(1.0, 2.0, 3.0)
	r_local = S.rotate_lab_to_local(geom, r_lab)
	r_back = S.rotate_local_to_lab(geom, r_local)
	@test r_back ≈ r_lab

	# Axis-aligned geometry: local z = lab z
	geom_z = S.Geometries.Cylindrical(; axis=SVector(0, 0, 1), radius=1.0, z=1.0 .. 5.0)
	R_z = S.rotation_local_to_lab(geom_z)
	@test R_z[:, 3] ≈ SVector(0, 0, 1)

	# Arbitrary axis roundtrip
	axis_arb = normalize(SVector(0.2, -0.7, 0.4))
	geom_arb = S.Geometries.Cylindrical(; axis=axis_arb, radius=0.3, z=1.0 .. 5.0)
	r_arb = SVector(1.1, -0.3, 2.7)
	r_local_arb = S.rotate_lab_to_local(geom_arb, r_arb)
	r_back_arb = S.rotate_local_to_lab(geom_arb, r_local_arb)
	@test r_back_arb ≈ r_arb
	@test r_local_arb[3] ≈ dot(r_arb, axis_arb)  # local z = axial component

	# Rotation matrices should match between Conical and Cylindrical for the same axis
	@testset for θ in [0, 0.1, π/2, deg2rad(175)]
		axis_y = normalize(SVector(sin(θ), 0, cos(θ)))
		geom_cyl = S.Geometries.Cylindrical(; axis=axis_y, radius=1.0, z=1.0 .. 5.0)
		geom_con = S.Geometries.Conical(; axis=axis_y, φj=0.1, z=1.0 .. 5.0)
		@test S.rotation_local_to_lab(geom_cyl) ≈ S.rotation_local_to_lab(geom_con)
	end
end

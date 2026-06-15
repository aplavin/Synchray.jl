@testitem "Parabolic geometry: natural_coords and is_inside" begin
	import Synchray as S

	# Canonical test geometry: R(z) = sqrt(z)
	# R(1)=1, R(4)=2, R(9)=3, R(25)=5, R(100)=10
	axis = SVector(0.0, 0.0, 1.0)
	geom = S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=0.5, z=1.0..100.0)

	@test S.geometry_axis(geom) ≈ axis

	@testset "natural_coords hardcoded" begin
		# On axis at z=4
		c = S.natural_coords(geom, S.FourPosition(0, 0, 0, 4))
		@test c.z ≈ 4.0
		@test c.ρ ≈ 0.0
		@test c.η ≈ 0.0

		# Off-axis at z=4, ρ=1 → η = 1/R(4) = 1/2
		c = S.natural_coords(geom, S.FourPosition(0, 1, 0, 4))
		@test c.z ≈ 4.0
		@test c.ρ ≈ 1.0
		@test c.η ≈ 0.5

		# Boundary at z=4, ρ=2 → η = 2/R(4) = 1
		c = S.natural_coords(geom, S.FourPosition(0, 2, 0, 4))
		@test c.z ≈ 4.0
		@test c.ρ ≈ 2.0
		@test c.η ≈ 1.0

		# Outside at z=4, ρ=3 → η = 3/2
		c = S.natural_coords(geom, S.FourPosition(0, 3, 0, 4))
		@test c.z ≈ 4.0
		@test c.ρ ≈ 3.0
		@test c.η ≈ 1.5

		# At z=25, ρ=3 → η = 3/R(25) = 3/5
		c = S.natural_coords(geom, S.FourPosition(0, 3, 0, 25))
		@test c.z ≈ 25.0
		@test c.ρ ≈ 3.0
		@test c.η ≈ 0.6

		# Val(:z) fast path
		@test S.natural_coords(geom, S.FourPosition(0, 1, 0, 4), Val(:z)) ≈ 4.0
	end

	@testset "is_inside hardcoded" begin
		@test S.is_inside(geom, S.FourPosition(0, 0, 0, 4)) == true    # on axis, in range
		@test S.is_inside(geom, S.FourPosition(0, 1, 0, 4)) == true    # ρ=1 < R(4)=2
		@test S.is_inside(geom, S.FourPosition(0, 2, 0, 4)) == true    # boundary
		@test S.is_inside(geom, S.FourPosition(0, 3, 0, 4)) == false   # ρ=3 > R(4)=2
		@test S.is_inside(geom, S.FourPosition(0, 0, 0, 0.5)) == false # z outside range
		@test S.is_inside(geom, S.FourPosition(0, 0, 0, 101)) == false # z outside range
	end

	@testset "η changes with z (unlike cylinder)" begin
		# Same ρ, different z → different η (because R(z) changes)
		c1 = S.natural_coords(geom, S.FourPosition(0, 1, 0, 4))
		c2 = S.natural_coords(geom, S.FourPosition(0, 1, 0, 25))
		@test c1.η ≈ 0.5   # 1/R(4) = 1/2
		@test c2.η ≈ 0.2   # 1/R(25) = 1/5
		@test c1.η != c2.η
	end
end


@testitem "Parabolic geometry: z_interval" begin
	import Synchray as S

	ray_at(x, y; ν=1.0) = S.RayZ(; x0=S.FourPosition(0.0, x, y, 0.0), k=ν, nz=16)

	# Helper for arbitrary-direction rays
	function arb_ray(r0, n̂; ν=1.0, nz=16)
		n̂ = normalize(n̂)
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		e1 = normalize(cross(up, n̂))
		e2 = cross(n̂, e1)
		S.Ray(S.FourPosition(0.0, r0...), S.photon_k(ν, n̂), e1, e2, nz)
	end

	# R(z) = sqrt(z): R(1)=1, R(4)=2, R(25)=5, R(100)=10
	geom = S.Geometries.Parabolic(; axis=SVector(0.0, 0.0, 1.0), R_ref=1.0, z_ref=1.0, a=0.5, z=1.0..100.0)

	@testset "axis = z, rays along z" begin
		# On axis: full truncation range
		@test S.z_interval(geom, ray_at(0.0, 0.0)) ≈ 1.0..100.0

		# ρ=1: inside for z ≥ 1 (R(1)=1=ρ), so entry at z=1
		seg = S.z_interval(geom, ray_at(1.0, 0.0))
		@test seg ≈ 1.0..100.0 rtol=1e-3

		# ρ=2: inside for z ≥ 4 (R(4)=2)
		seg = S.z_interval(geom, ray_at(2.0, 0.0))
		@test seg ≈ 4.0..100.0 rtol=1e-3

		# ρ=5: inside for z ≥ 25
		seg = S.z_interval(geom, ray_at(5.0, 0.0))
		@test seg ≈ 25.0..100.0 rtol=1e-3

		# ρ=11: would need z ≥ 121 > 100 → miss
		@test S.z_interval(geom, ray_at(11.0, 0.0)) |> isempty
	end

	@testset "perpendicular ray" begin
		# Ray along x at z=4: R(4)=2, inside when |s| ≤ 2
		seg = S.z_interval(geom, arb_ray(SA[0, 0, 4], SA[1, 0, 0]))
		@test seg ≈ -2.0..2.0 rtol=1e-3

		# Ray along x at z=25: R(25)=5, inside when |s| ≤ 5
		seg = S.z_interval(geom, arb_ray(SA[0, 0, 25], SA[1, 0, 0]))
		@test seg ≈ -5.0..5.0 rtol=1e-3

		# Ray along x at z=0.5 (outside truncation): miss
		@test S.z_interval(geom, arb_ray(SA[0, 0, 0.5], SA[1, 0, 0])) |> isempty
	end

	@testset "tangent ray (single contact point)" begin
		# 45° ray from origin: z(s) = s/√2, ρ(s) = s/√2
		# Inside when s/√2 ≤ √(s/√2), i.e. s ≤ √2
		# But z ≥ 1 requires s ≥ √2
		# Contact at exactly one point → effectively empty
		seg = S.z_interval(geom, arb_ray(SA[0, 0, 0], SA[0, 1/√2, 1/√2]))
		@test S.width(seg) < 0.2  # effectively empty or very narrow (conservative bisection widens slightly)
	end

	@testset "different shape exponents" begin
		# a=1 (conical): R(z) = z * (R_ref/z_ref)
		geom_cone = S.Geometries.Parabolic(; axis=SVector(0.0, 0.0, 1.0), R_ref=0.5, z_ref=1.0, a=1.0, z=1.0..10.0)
		# On axis
		@test S.z_interval(geom_cone, ray_at(0.0, 0.0)) ≈ 1.0..10.0
		# ρ=2: inside when 2 ≤ 0.5z → z ≥ 4
		seg = S.z_interval(geom_cone, ray_at(2.0, 0.0))
		@test seg ≈ 4.0..10.0 rtol=1e-3

		# a=0.25: R(z) = z^0.25 (slow widening)
		geom_slow = S.Geometries.Parabolic(; axis=SVector(0.0, 0.0, 1.0), R_ref=1.0, z_ref=1.0, a=0.25, z=1.0..100.0)
		# ρ=2: inside when 2 ≤ z^0.25 → z ≥ 16
		seg = S.z_interval(geom_slow, ray_at(2.0, 0.0))
		@test seg ≈ 16.0..100.0 rtol=1e-3
	end

	@testset "tilted axis" begin
		θ = 0.2
		axis_tilted = SVector(sin(θ), 0.0, cos(θ))
		geom_tilted = S.Geometries.Parabolic(; axis=axis_tilted, R_ref=1.0, z_ref=1.0, a=0.5, z=1.0..100.0)

		# On-axis ray should still see the full truncation range
		seg = S.z_interval(geom_tilted, ray_at(0.0, 0.0))
		# The z-axis projects onto the tilted axis with cos(θ), so intersection covers
		# z_lab ∈ [1/cos(θ), 100/cos(θ)], but only if the z-axis is inside the envelope.
		# For a parabolic envelope, the z-axis is inside when ρ_to_axis < R(z_along_axis).
		# ρ = z_lab * sin(θ), z_along = z_lab * cos(θ)
		# Inside when z_lab * sin(θ) ≤ (z_lab * cos(θ))^0.5
		# → z_lab * sin(θ) ≤ z_lab^0.5 * cos(θ)^0.5
		# → z_lab^0.5 ≤ cos(θ)^0.5 / sin(θ)
		# → z_lab ≤ cos(θ) / sin(θ)^2
		z_exit = cos(θ) / sin(θ)^2
		@test seg ≈ 1.02034 .. z_exit / cos(θ) rtol=0.02
	end
end


@testitem "Parabolic geometry: consistency with Conical at a=1" begin
	import Synchray as S

	ray_at(x, y; ν=1.0) = S.RayZ(; x0=S.FourPosition(0.0, x, y, 0.0), k=ν, nz=16)

	φj = 0.05
	z = 1e-3 .. 5.0
	axis = SVector(0.0, 0.0, 1.0)
	z_ref = 1.0

	geom_cone = S.Geometries.Conical(; axis, φj, z)
	geom_para = S.Geometries.Parabolic(; axis, R_ref=z_ref * tan(φj), z_ref, a=1.0, z)

	@testset "natural_coords match" begin
		for x4 in [S.FourPosition(0, 0, 0, 1), S.FourPosition(0, 0.01, 0, 2), S.FourPosition(0, 0.03, 0.01, 3)]
			cc = S.natural_coords(geom_cone, x4)
			cp = S.natural_coords(geom_para, x4)
			@test cc.z ≈ cp.z
			@test cc.ρ ≈ cp.ρ
			@test cc.η ≈ cp.η rtol=1e-10
		end
	end

	@testset "is_inside match" begin
		for x4 in [
			S.FourPosition(0, 0, 0, 1),       # on axis
			S.FourPosition(0, 0.01, 0, 2),     # inside
			S.FourPosition(0, 0.3, 0, 2),      # outside (ρ > z*tan(φ))
			S.FourPosition(0, 0, 0, 6),        # z outside range
		]
			@test S.is_inside(geom_cone, x4) == S.is_inside(geom_para, x4)
		end
	end

	@testset "z_interval match" begin
		for (x, y) in [(0.0, 0.0), (0.05, 0.0), (0.1, 0.0), (0.5, 0.0)]
			seg_c = S.z_interval(geom_cone, ray_at(x, y))
			seg_p = S.z_interval(geom_para, ray_at(x, y))
			if isempty(seg_c)
				@test isempty(seg_p)
			else
				@test leftendpoint(seg_p) ≈ leftendpoint(seg_c) rtol=0.01
				@test rightendpoint(seg_p) ≈ rightendpoint(seg_c) rtol=0.01
			end
		end
	end

	@testset "z_interval match with tilted axis" begin
		θ = 0.2
		axis_t = SVector(sin(θ), 0.0, cos(θ))
		geom_cone_t = S.Geometries.Conical(; axis=axis_t, φj, z)
		geom_para_t = S.Geometries.Parabolic(; axis=axis_t, R_ref=z_ref * tan(φj), z_ref, a=1.0, z)

		for (x, y) in [(0.0, 0.0), (0.05, 0.0), (0.0, 0.05), (0.5, 0.0)]
			seg_c = S.z_interval(geom_cone_t, ray_at(x, y))
			seg_p = S.z_interval(geom_para_t, ray_at(x, y))
			if isempty(seg_c)
				@test isempty(seg_p)
			else
				@test leftendpoint(seg_p) ≈ leftendpoint(seg_c) rtol=0.01
				@test rightendpoint(seg_p) ≈ rightendpoint(seg_c) rtol=0.01
			end
		end
	end
end


@testitem "Parabolic geometry: constructor validation" begin
	import Synchray as S

	axis = SVector(0.0, 0.0, 1.0)

	# a > 1 should fail
	@test_throws AssertionError S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=1.5, z=1.0..10.0)

	# a = 0 should fail
	@test_throws AssertionError S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=0.0, z=1.0..10.0)

	# a = -1 should fail
	@test_throws AssertionError S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=-1.0, z=1.0..10.0)

	# a = 1 should work
	geom = S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=1.0, z=1.0..10.0)
	@test geom.a == 1.0

	# a = 0.5 should work
	geom = S.Geometries.Parabolic(; axis, R_ref=1.0, z_ref=1.0, a=0.5, z=1.0..10.0)
	@test geom.a == 0.5
end

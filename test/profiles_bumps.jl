@testitem "Inverse bump patterns" begin
	import Synchray as S

	@testset "TophatBump inverse - exact pointwise reciprocal" begin
		b = S.Patterns.TophatBump(2.0)
		b_inv = inv(b)

		# Exact pointwise inverse at all χ
		for χ in [0.0, 0.5, 0.9, 1.0, 1.5, 5.0]
			@test b(χ) * b_inv(χ) ≈ 1.0
		end

		# Specific values
		@test b(0.5) ≈ 2.0  # Inside: χ < 1
		@test b_inv(0.5) ≈ 0.5  # Exact reciprocal
		@test b(1.5) ≈ 1.0  # Outside: χ ≥ 1
		@test b_inv(1.5) ≈ 1.0  # Reciprocal of 1 is 1

		# Double inversion recovers values
		@test inv(inv(b))(0.5) ≈ b(0.5)
		@test inv(inv(b))(1.5) ≈ b(1.5)
	end

	@testset "GaussianBump inverse - exact pointwise reciprocal" begin
		g = S.Patterns.GaussianBump(f_peak=3.0, χ_threshold=16.0)
		g_inv = inv(g)

		# Exact pointwise inverse at all χ
		for χ in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 20.0]
			@test g(χ) * g_inv(χ) ≈ 1.0  atol=1e-15
		end

		# Specific values
		@test g(0.0) ≈ 3.0
		@test g_inv(0.0) ≈ 1/3
		@test g(20.0) ≈ 1.0  # Beyond threshold
		@test g_inv(20.0) ≈ 1.0

		# Double inversion recovers values
		for χ in [0.0, 1.0, 5.0, 20.0]
			@test inv(inv(g))(χ) ≈ g(χ)
		end
	end

	@testset "Identity case" begin
		# f_peak = 1.0 means no enhancement, so inverse should also be identity
		b_id = S.Patterns.TophatBump(1.0)
		b_inv = inv(b_id)
		for χ in [0.0, 0.5, 1.5]
			@test b_id(χ) ≈ 1.0
			@test b_inv(χ) ≈ 1.0
		end

		g_id = S.Patterns.GaussianBump(f_peak=1.0)
		g_inv = inv(g_id)
		for χ in [0.0, 1.0, 5.0]
			@test g_id(χ) ≈ 1.0
			@test g_inv(χ) ≈ 1.0
		end
	end

	@testset "Type preservation" begin
		# Float32
		b32 = S.Patterns.TophatBump(2.0f0)
		b_inv32 = inv(b32)
		@test b_inv32(0.5f0) isa Float32
		@test b32(0.5f0) * b_inv32(0.5f0) ≈ 1.0f0

		g32 = S.Patterns.GaussianBump(f_peak=2.0f0, χ_threshold=16.0f0)
		g_inv32 = inv(g32)
		@test g_inv32(1.0f0) isa Float32
		@test g32(1.0f0) * g_inv32(1.0f0) ≈ 1.0f0
	end

	@testset "Suppression profiles (f_peak < 1)" begin
		# Test that suppressions work correctly
		g_suppress = S.Patterns.GaussianBump(f_peak=0.5)  # Suppression
		g_enhance = inv(g_suppress)  # Inverted suppression = enhancement

		@test g_suppress(0.0) ≈ 0.5  # Suppresses to 50% at center
		@test g_enhance(0.0) ≈ 2.0   # Enhances by 2x at center
		@test g_suppress(0.0) * g_enhance(0.0) ≈ 1.0  # Exact inverse
	end

	@testset "Integration with knot patterns" begin
		# Verify inverted profiles work when used as knot profiles
		knot_enhance = S.Patterns.EllipsoidalKnot(
			x_c0 = S.FourPosition(0.0, 0.0, 0.0, 2.0),
			u = S.FourVelocity(SVector(0.0, 0.0, 0.1)),
			sizing = S.Patterns.FixedSizing(0.2, 0.5),
			profile = S.Patterns.GaussianBump(3.0),
		)

		knot_suppress = S.Patterns.EllipsoidalKnot(
			x_c0 = knot_enhance.x_c0,
			u = knot_enhance.u,
			sizing = knot_enhance.sizing,
			profile = inv(knot_enhance.profile),  # Inverted profile
		)

		# Create a simple geometry for testing
		geom = S.Geometries.Conical(; axis = SVector(0, 0, 1), φj = 0.05, z = 1e-3..10)

		# At the center, enhancement and suppression should be exact inverses
		x4c = knot_enhance.x_c0
		factor_enhance = knot_enhance(geom, x4c, 1.0)
		factor_suppress = knot_suppress(geom, x4c, 1.0)

		@test factor_enhance ≈ 3.0
		@test factor_suppress ≈ 1/3
		@test factor_enhance * factor_suppress ≈ 1.0
	end
end

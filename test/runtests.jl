using TestItems
using TestItemRunner
@run_package_tests


@testitem "Minkowski + Doppler conventions" begin
	import Synchray as S

	@testset "arithmetics" begin
		a = S.FourPosition(1, 2, 3, 4)
		b = S.FourPosition(0.5, -1.0, 2.0, 3.0)
		@test a + b === S.FourPosition(1.5, 1.0, 5.0, 7.0)
		@test a - b === S.FourPosition(0.5, 3.0, 1.0, 1.0)
		@test 2a === S.FourPosition(2, 4, 6, 8)
	end

	@testset "minkowski_dot basics" begin
		a = S.FourPosition(2.0, 1.0, 3.0, -4.0)
		b = S.FourPosition(-1.0, 2.0, 0.5, 7.0)
		@test S.minkowski_dot(a, b) ≈ S.minkowski_dot(b, a)
		@test S.minkowski_dot(a, a) ≈ -(a.t^2) + a.x^2 + a.y^2 + a.z^2
	end

	@testset "FourVelocity normalization" begin
		β = SVector(0.3, -0.4, 0.1)
		u = S.FourVelocity(β)
		@test S.minkowski_dot(u, u) ≈ -1 atol=1e-12
		@test u.t ≈ inv(sqrt(1 - dot(β, β)))
	end

	@testset "photon_k is null" begin
		ν = 2.5
		n = normalize(SVector(0.2, -0.3, 0.7))
		k = S.photon_k(ν, n)
		@test S.minkowski_dot(k, k) ≈ 0 atol=1e-12
	end

	@testset "Doppler boosting δ convention" begin
		# Convention locked in by this test:
		# doppler_factor(u, n) == δ = ν_obs / ν_comoving
		# For fast motion toward the observer along the photon direction (β ⋅ n > 0): δ >> 1.
		β = 0.99
		u_toward = S.FourVelocity(SVector(0.0, 0.0, β))
		n = SVector(0.0, 0.0, 1.0)

		δ_toward = S.doppler_factor(u_toward, n)
		@test δ_toward > 1
		@test δ_toward > 10

		u_away = S.FourVelocity(SVector(0.0, 0.0, -β))
		δ_away = S.doppler_factor(u_away, n)
		@test δ_away < 1
	end
end


@testitem "Uniform slab transfer" begin
    import Synchray as S

	L = 3
	j0 = 0.7
	a0 = 1.3

	ray = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2, nz=2_000)

	cases = [
		(; u0=S.FourVelocity(SVector(0, 0, 0))),
		(; u0=S.FourVelocity(SVector(0, 0, 0.3))),
		(; u0=S.FourVelocity(SVector(0.3, 0, 0.5))),
	]

	@testset for (;u0) in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset "limits" begin
			@testset "optically thin" begin
				# For τ = α D L ≪ 1: I′ ≈ j L′ and I = I′/D^3 ⇒ I ≈ j L / D^2
				αthin = 1e-6
				slab = S.UniformSlab(0..L, u0, j0, αthin)
				Iν_num = S.render(ray, slab)
				Iν_exact = j0 * L * (δ^2)
				@test Iν_num ≈ Iν_exact rtol=2e-3
				@test S.render(ray, slab, S.OpticalDepth()) ≈ 2e-6 rtol=0.5
			end

			@testset "very optically thick (α ≫ j)" begin
				# For τ ≫ 1: I′ → S = j/α and I = S/D^3.
				αthick = 1e5
				slab = S.UniformSlab(0..L, u0, j0, αthick)
				Iν_num = S.render(ray, slab)
				Iν_exact = (j0 / αthick) * (δ^3)
				@test Iν_num ≈ Iν_exact rtol=2e-3
				@test S.render(ray, slab, S.OpticalDepth()) ≈ 2e5 rtol=0.5
			end
		end

		slab = S.UniformSlab(0.0..L, u0, j0, a0)
		Iν_num = S.render(ray, slab)
		τ_num = S.render(ray, slab, S.OpticalDepth())
		α_num = S.render(ray, slab, S.SpectralIndex())

		Iν_τ_num = S.render(ray, slab, (S.Intensity(), S.OpticalDepth()))
		@test Iν_τ_num[1] ≈ Iν_num
		@test Iν_τ_num[2] ≈ τ_num

		Iν_α_num = S.render(ray, slab, (S.Intensity(), S.SpectralIndex()))
		@test Iν_α_num[1] ≈ Iν_num
		@test Iν_α_num[2] ≈ α_num

		# Analytic solution as a cross-check:
		# - Solve transfer in slab comoving frame: dI'/ds' = j0 - a0 I'
		# - Convert thickness to comoving path length via D = ν'/ν = 1/δ
		#   for a ray along +z (n = ẑ): δ = doppler_factor(u, ẑ)
		# - Transform specific intensity back: I_ν = δ^{3} I'_{ν'}
		L′ = L / δ
		I′ = (j0 / a0) * (1 - exp(-a0 * L′))
		Iν_exact = I′ * (δ^3)

		@test Iν_num ≈ Iν_exact rtol=2e-3
		@test τ_num ≈ a0 * L′ rtol=2e-3

		# For frequency-independent comoving (j0, a0), the observed I_ν is ν-independent ⇒ α = 0.
		@test α_num ≈ 0 atol=2e-6
	end
end


@testitem "Uniform sphere flux (thin to thick)" begin
	import Synchray as S
    using RectiGrids
	using Accessors: @set

	R = 1.3
	j0 = 0.7
	center = S.FourPosition(0, 0, 0, 0)
	αs = (0, 1.2, 12)

	ray = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2, nz=512)

	cases = (
		S.FourVelocity(SVector(0, 0, 0)),
		S.FourVelocity(SVector(0, 0, 0.3)),
		S.FourVelocity(SVector(0.5, 0, 0.3)),
	)

	I_exact(j0, α0, δ, ℓ) = α0 == 0 ? (j0 * ℓ) * (δ^2) : (j0 / α0) * (1 - exp(-α0 * ℓ / δ)) * (δ^3)
	F_exact(j0, α0, δ, R) = α0 == 0 ? (4π / 3) * j0 * R^3 * (δ^2) : begin
		S0 = j0 / α0
		τ0 = 2 * α0 * R / δ
		(π * R^2 * S0 * (δ^3)) * (1 - (2 / (τ0^2)) * (1 - (1 + τ0) * exp(-τ0)))
	end

	@testset for u0 in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset for α0 in αs
			sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=α0)

			# Per-ray analytic check
			for b in (0, 0.4R, 0.8R)
				ℓ = 2 * sqrt(R^2 - b^2)
				I_num = S.render((@set ray.x0.x = b), sphere)
				@test I_num ≈ I_exact(j0, α0, δ, ℓ) rtol=6e-3
			end

			# Missed rays should return zero/nan
			mray = @set ray.x0.x = 2R
			@test S.render(mray, sphere) == 0
			@test S.render(mray, sphere, S.OpticalDepth()) == 0
			@test S.render(mray, sphere, S.SpectralIndex()) |> isnan
			@test isequal(S.render(mray, sphere, (S.Intensity(), S.SpectralIndex())), (0., NaN))

			# Total (image-plane) flux from a 2D grid: ∫ I(x,y) dx dy
			xs = range(-R..R, 151)
			cam = S.CameraZ(; xys=grid(SVector, xs, xs), ray.nz, ν=2, t=0)
			img = S.render(cam, sphere)
			dx = step(xs)
			F_num = sum(img) * dx^2
			@test F_num ≈ F_exact(j0, α0, δ, R) rtol=(α0 == 0 ? 1.6e-2 : 1.2e-2)
		end
	end
end


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

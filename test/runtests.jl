using TestItems
using TestItemRunner
@run_package_tests


@testitem "Uniform slab transfer" begin
    import Synchray as S

	νcam = 2.0
	L = 3.0
	j0 = 0.7
	a0 = 1.3

	cases = [
		(; u0=S.FourVelocity(SVector(0.0, 0.0, 0.0))),
		(; u0=S.FourVelocity(SVector(0.0, 0.0, 0.3))),
		(; u0=S.FourVelocity(SVector(0.3, 0.0, 0.5))),
	]

	@testset for (;u0) in cases
		D = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset "limits" begin
			@testset "optically thin" begin
				# For τ = α D L ≪ 1: I′ ≈ j L′ and I = I′/D^3 ⇒ I ≈ j L / D^2
				αthin = 1e-6
				slab = S.UniformSlab(0.0..L, u0, j0, αthin)
				Iν_num = S.integrate_ray(slab, SVector(0, 0); νcam, t0=0.0, nz=20_000)
				Iν_exact = j0 * L / (D^2)
				@test Iν_num ≈ Iν_exact rtol=2e-3
			end

			@testset "very optically thick (α ≫ j)" begin
				# For τ ≫ 1: I′ → S = j/α and I = S/D^3.
				αthick = 1e5
				slab = S.UniformSlab(0.0..L, u0, j0, αthick)
				Iν_num = S.integrate_ray(slab, SVector(0, 0); νcam, t0=0.0, nz=20_000)
				Iν_exact = (j0 / αthick) / (D^3)
				@test Iν_num ≈ Iν_exact rtol=2e-3
			end
		end

		slab = S.UniformSlab(0.0..L, u0, j0, a0)
		Iν_num = S.integrate_ray(slab, SVector(0, 0); νcam, t0=0.0, nz=20_000)

		# Analytic solution as a cross-check:
		# - Solve transfer in slab comoving frame: dI'/ds' = j0 - a0 I'
		# - Convert thickness to comoving path length via Doppler factor D = γ(1-β·n)
		#   for a ray along +z (n = ẑ): D = doppler_factor(u, ẑ)
		# - Transform specific intensity back: I_ν = D^{-3} I'_{ν'}
		L′ = D * L
		I′ = (j0 / a0) * (1 - exp(-a0 * L′))
		Iν_exact = I′ / (D^3)

		@test Iν_num ≈ Iν_exact rtol=2e-3
	end
end


@testitem "Uniform sphere flux (thin to thick)" begin
	import Synchray as S
    using RectiGrids

	R = 1.3
	j0 = 0.7
	νcam = 2.0
	center = S.FourPosition(0.0, 0.0, 0.0, 0.0)
	αs = (0.0, 1.2, 12.0)

	cases = (
		S.FourVelocity(SVector(0.0, 0.0, 0.0)),
		S.FourVelocity(SVector(0.0, 0.0, 0.3)),
		S.FourVelocity(SVector(0.5, 0.0, 0.3)),
	)

	I_exact(j0, α0, D, ℓ) = α0 == 0 ? (j0 * ℓ) / (D^2) : (j0 / α0) * (1 - exp(-α0 * D * ℓ)) / (D^3)
	F_exact(j0, α0, D, R) = α0 == 0 ? (4π / 3) * j0 * R^3 / (D^2) : begin
		S0 = j0 / α0
		τ0 = 2 * α0 * D * R
		(π * R^2 * S0 / (D^3)) * (1 - (2 / (τ0^2)) * (1 - (1 + τ0) * exp(-τ0)))
	end

	@testset for u0 in cases
		D = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset for α0 in αs
			sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=α0)

			# Per-ray analytic check
			for b in (0.0, 0.4R, 0.8R)
				ℓ = 2 * sqrt(R^2 - b^2)
				I_num = S.integrate_ray(sphere, SVector(b, 0.0); νcam, t0=0.0, nz=256)
				@test I_num ≈ I_exact(j0, α0, D, ℓ) rtol=6e-3
			end

			# Missed rays should return zero
			@test S.integrate_ray(sphere, SVector(2R, 0.0); νcam, t0=0.0, nz=64) == 0

			# Total (image-plane) flux from a 2D grid: ∫ I(x,y) dx dy
			xs = range(-R..R, 151)
			cam = S.OrthoCamera(; xys=grid(SVector, xs, xs), nz=256, ν=νcam, t=0.0)
			img = S.render(cam, sphere)
			dx = step(xs)
			F_num = sum(img) * dx^2
			@test F_num ≈ F_exact(j0, α0, D, R) rtol=(α0 == 0 ? 1.6e-2 : 1.2e-2)
		end
	end
end


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

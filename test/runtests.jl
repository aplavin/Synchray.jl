using TestItems
using TestItemRunner
@run_package_tests


@testitem "Invariant slab transfer" begin
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
		slab = S.UniformSlab(0.0..L, u0, j0, a0)
		Iν_num = S.integrate_ray(slab, SVector(0, 0); νcam, t0=0.0, nz=20_000)

		# Analytic solution as a cross-check:
		# - Solve transfer in slab comoving frame: dI'/ds' = j0 - a0 I'
		# - Convert thickness to comoving path length via Doppler factor D = γ(1-β·n)
		#   for a ray along +z (n = ẑ): D = u.t - u.z
		# - Transform specific intensity back: I_ν = D^{-3} I'_{ν'}
		D = u0.t - u0.z
		L′ = D * L
		I′ = (j0 / a0) * (1 - exp(-a0 * L′))
		Iν_exact = I′ / (D^3)

		@test Iν_num ≈ Iν_exact rtol=2e-3
	end
end


@testitem "Optically thin uniform sphere flux" begin
	import Synchray as S
    using RectiGrids

	R = 1.3
	j0 = 0.7
	νcam = 2.0
	center = S.FourPosition(0.0, 0.0, 0.0, 0.0)
	α0 = 0.0

	cases = [
		(; u0=S.FourVelocity(SVector(0.0, 0.0, 0.0))),
		(; u0=S.FourVelocity(SVector(0.0, 0.0, 0.3))),
		(; u0=S.FourVelocity(SVector(0.5, 0.0, 0.3))),
	]

	@testset for (; u0) in cases
		sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=α0)
		D = u0.t - u0.z

		# Per-ray analytic check: I = j0 * ℓ / D^2, ℓ = 2√(R^2 - b^2)
		for b in (0.0, 0.4R, 0.8R)
			ℓ = 2 * sqrt(R^2 - b^2)
			I_exact = j0 * ℓ / (D^2)
			I_num = S.integrate_ray(sphere, SVector(b, 0.0); νcam, t0=0.0, nz=4_000)
			@test I_num ≈ I_exact rtol=4e-3
		end

		# Missed rays should return zero
		@test S.integrate_ray(sphere, SVector(2R, 0.0); νcam, t0=0.0, nz=64) == 0

		# Total (image-plane) flux from a 2D grid: ∫ I(x,y) dx dy
		# For uniform optically thin sphere: (4π/3) * j0 * R^3 / D^2.
		xs = range(-R..R, 101)
		cam = S.OrthoCamera(xys=grid(SVector, xs, xs), nz=2_000, ν=νcam, t=0.0)
		img = S.render(cam, sphere)
		dx = step(xs)
		F_num = sum(img) * dx * dx
		F_exact = (4π / 3) * j0 * R^3 / (D^2)
		@test F_num ≈ F_exact rtol=1.2e-2
	end
end


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

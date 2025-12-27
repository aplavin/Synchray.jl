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


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

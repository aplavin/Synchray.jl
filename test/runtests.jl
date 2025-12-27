using TestItems
using TestItemRunner
@run_package_tests


@testitem "Invariant slab transfer" begin
    import Synchray as S

	νcam = 2.0
	L = 3.0
	j0 = 0.7
	a0 = 1.3
	slab = S.ConstantSlab(0.0..L, S.FourVelocity(SVector(0,0,0)), j0, a0)

	Iν_num = S.integrate_ray(slab, 0.0, 0.0; νcam, t0=0.0, nz=20_000)
	Δλ = L / S.photon_k(νcam).z
	Jinv = j0 / (νcam^2)
	Ainv = a0 * νcam
	Iinv_exact = (Jinv / Ainv) * (1 - exp(-Ainv * Δλ))
	Iν_exact = νcam^3 * Iinv_exact

	@test isapprox(Iν_num, Iν_exact; rtol=2e-3, atol=0)
end


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

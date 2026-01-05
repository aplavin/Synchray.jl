@testitem "Polarized synchrotron coeffs (intrinsic)" begin
	import Synchray as S
	using Test
	using Accessors

	ν = 2.0
	x4 = S.FourPosition(0.0, 0.0, 0.0, 0.0)
	k′ = S.photon_k(ν, SVector(0.0, 0.0, 1.0))

	u0 = S.FourVelocity(SVector(0.0, 0.0, 0.0))
	ne0 = 1.3
	B0 = 0.9
	z = 0.0..1.0

	cases = (
		S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.3),
		S.AnisotropicPowerLawElectrons(; p=2.5, η=0.5, Cj=0.7, Ca=0.3),
	)

	@testset for model in cases
		slab = S.UniformSynchrotronSlab(; z, u0, ne0, B0=SVector(B0, 0.0, 0.0), model)
		(j, a) = S.emissivity_absorption_polarized(slab, x4, k′)
		(jI, αI) = S.emissivity_absorption(slab, x4, k′)

		@test isfinite(j.perp) && isfinite(j.par)
		@test isfinite(a.perp) && isfinite(a.par)
		@test j.perp ≥ 0 && j.par ≥ 0
		@test a.perp ≥ 0 && a.par ≥ 0
		@test (j.perp + j.par) ≈ jI rtol=5e-15 atol=0
		@test (a.perp + a.par) ≈ αI rtol=5e-15 atol=0
	end

	@testset "FullyTangled depolarizes" begin
		model = S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.3)
		slab = S.UniformSynchrotronSlab(; z, u0, ne0, B0=S.FullyTangled(B0), model)
		(j, a) = S.emissivity_absorption_polarized(slab, x4, k′)
		@test j.perp ≈ j.par
		@test a.perp ≈ a.par
		(jI, αI) = S.emissivity_absorption(slab, x4, k′)
		@test (j.perp + j.par) ≈ jI rtol=5e-15 atol=0
		@test (a.perp + a.par) ≈ αI rtol=5e-15 atol=0
	end

	@testset "Mode↔Stokes helper sanity" begin
		m = S.ModePerpPar(3.0, 1.0)
		IQ = S.stokes_IQ(m)
		@test IQ.I ≈ 4.0
		@test IQ.Q ≈ 2.0

		m2 = S.modes_from_IQ(IQ.I, IQ.Q)
		@test m2.perp ≈ m.perp
		@test m2.par ≈ m.par
	end

	@testset "Polarized invariants" begin
		model = S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.3)
		slab = S.UniformSynchrotronSlab(; z, u0, ne0, B0=SVector(B0, 0.0, 0.0), model)
		(j, a) = S.emissivity_absorption_polarized(slab, x4, k′)
		(Jinv, Ainv) = S.emissivity_absorption_polarized_invariant(slab, x4, k′)
		@test Jinv ≈ j / (ν^2)
		@test Ainv ≈ a * ν
	end
end

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
		@test (a.perp + a.par)/2 ≈ αI rtol=5e-15 atol=0
	end

	@testset "FullyTangled depolarizes" begin
		model = S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.3)
		slab = S.UniformSynchrotronSlab(; z, u0, ne0, B0=S.FullyTangled(B0), model)
		(j, a) = S.emissivity_absorption_polarized(slab, x4, k′)
		@test j.perp ≈ j.par
		@test a.perp ≈ a.par
		(jI, αI) = S.emissivity_absorption(slab, x4, k′)
		@test (j.perp + j.par) ≈ jI rtol=5e-15 atol=0
		@test (a.perp + a.par)/2 ≈ αI rtol=5e-15 atol=0
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

@testitem "Stokes QU rotation" begin
	import Synchray as S
	using Test

	@testset "R then R^{-1}" begin
		QU0 = SVector(0.7, -1.2)
		χ = 0.37
		R = S.stokes_QU_rotation(χ)
		Rinv = S.stokes_QU_rotation(-χ)
		QU1 = S.rotate_QU(R, QU0)
		QU2 = S.rotate_QU(Rinv, QU1)
		@test QU2 ≈ QU0
	end

	@testset "45° rotation sign check" begin
		# For χ = 45° (π/4), the Stokes rotation is by 2χ = 90°:
		# Q' = U,  U' = -Q.
		χ = π / 4
		R = S.stokes_QU_rotation(χ)
		@test R ≈ (@SMatrix [0 1; -1 0])

		QU0 = SVector(3, 5)
		QU1 = S.rotate_QU(R, QU0)
		@test QU1[1] ≈ QU0[2]
		@test QU1[2] ≈ -QU0[1]
	end

	@testset "Simple geometry gives identity" begin
		ν = 2
		k = S.photon_k(ν, SVector(0, 0, 1))
		u = S.FourVelocity(SVector(0, 0, 0))
		(n′, e1′, e2′) = S.comoving_screen_basis(u, k)

		(e_par, e_perp) = S.linear_polarization_basis_from_B(n′, SVector(1, 0, 0))
		R = S.stokes_QU_rotation(e1′, e2′, e_par)
		@test R ≈ SMatrix{2,2}(I)
	end

	@testset "Known nonzero rotation in screen plane" begin
		ν = 2
		k = S.photon_k(ν, SVector(0, 0, 1))
		u = S.FourVelocity(SVector(0, 0, 0))
		(n′, e1′, e2′) = S.comoving_screen_basis(u, k)

		# Choose B′ in the comoving screen plane at a known angle χ from e1′.
		χ = 0.41
		B′ = cos(χ) * e1′ + sin(χ) * e2′
		(e_par, _) = S.linear_polarization_basis_from_B(n′, B′)
		@test dot(e_par, B′) ≈ norm(B′)

		Rgeom = S.stokes_QU_rotation(e1′, e2′, e_par)
		Rana = S.stokes_QU_rotation(χ)
		@test Rgeom ≈ Rana
	end

	@testset "Boosted screen basis + identity when B′ projects to e1′" begin
		ν = 2
		k = S.photon_k(ν, SVector(0, 0, 1))
		u = S.FourVelocity(SVector(0.3, 0.2, 0.5))
		(n′, e1′, e2′) = S.comoving_screen_basis(u, k)

		# Basis should remain orthonormal and right-handed in the comoving frame.
		@test norm(n′) ≈ 1
		@test norm(e1′) ≈ 1
		@test norm(e2′) ≈ 1
		@test dot(n′, e1′) ≈ 0  atol=√eps(1.)
		@test dot(n′, e2′) ≈ 0  atol=√eps(1.)
		@test dot(e1′, e2′) ≈ 0  atol=√eps(1.)
		@test cross(n′, e1′) ≈ e2′

		# Construct B′ with a nonzero component along n′; projection should still align with e1′.
		B′ = e1′ + 0.2 * n′
		(e_par, _) = S.linear_polarization_basis_from_B(n′, B′)
		R = S.stokes_QU_rotation(e1′, e2′, e_par)
		@test R ≈ SMatrix{2,2}(I)
	end

	@testset "Edge cases (expected failures / non-finite outputs)" begin
		ν = 2
		u0 = S.FourVelocity(SVector(0, 0, 0))

		@testset "Ray parallel to camera axes" begin
			kx = S.photon_k(ν, SVector(1, 0, 0))
			ky = S.photon_k(ν, SVector(0, 1, 0))
			@test_throws AssertionError S.comoving_screen_basis(u0, kx)
			@test_throws AssertionError S.comoving_screen_basis(u0, ky)
		end

		@testset "Non-unit / invalid ray directions" begin
			kbad = S.photon_k(ν, SVector(2, 0, 0))
			kzero = S.photon_k(ν, SVector(0, 0, 0))
			@test_throws AssertionError S.comoving_screen_basis(u0, kbad)
			@test_throws AssertionError S.comoving_screen_basis(u0, kzero)
		end

		@testset "Magnetic field degeneracies" begin
			n = SVector(0, 0, 1)
			@test_throws AssertionError S.linear_polarization_basis_from_B(n, SVector(0, 0, 0))

			(e_par, e_perp) = S.linear_polarization_basis_from_B(n, n)
			@test norm(e_par) ≈ 1
			@test norm(e_perp) ≈ 1
		end

		@testset "Arbitrary basis for zero direction" begin
			(e1, e2) = S._arbitrary_screen_basis(SVector(0, 0, 0))
			@test any(!isfinite, e1) || any(!isfinite, e2)
		end
	end
end

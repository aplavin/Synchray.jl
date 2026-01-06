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


@testitem "Linear-polarized uniform slab transfer" begin
	import Synchray as S
	using Test

	struct ConstPolarizedSlab{TZ,TU,TB,TJ,TA} <: S.AbstractSynchrotronMedium
		z::TZ
		u0::TU
		B0::TB
		j_mode::TJ
		a_mode::TA
	end

	S.z_interval(obj::ConstPolarizedSlab, ray::S.RayZ) = obj.z
	S.four_velocity(obj::ConstPolarizedSlab, x4) = obj.u0
	S.magnetic_field(obj::ConstPolarizedSlab, x4) = obj.B0
	S.emissivity_absorption_polarized(obj::ConstPolarizedSlab, x4, k′) = (obj.j_mode, obj.a_mode)
	S.emissivity_absorption(obj::ConstPolarizedSlab, x4, k′) = (obj.j_mode.perp + obj.j_mode.par, (obj.a_mode.perp + obj.a_mode.par) / 2)

	ν = 2.0
	# Pick χ so that e_perp (the +Q axis) is at angle χ from the camera x-axis.
	χ = 0.31
	B′ = SVector(-sin(χ), cos(χ), 0.0)
	
	z = 0.0 .. 1.0
	u0 = S.FourVelocity(SVector(0.0, 0.0, 0.0))

	# Choose coefficients so Δτ is safely above the linear threshold.
	j_mode = S.ModePerpPar(1.7, 0.4)
	a_mode = S.ModePerpPar(2.0, 1.0)
	obj = ConstPolarizedSlab(z, u0, B′, j_mode, a_mode)

	ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=S.photon_k(ν, SVector(0.0, 0.0, 1.0)), nz=40)

	S_cam = S.render(ray, obj, S.IntensityIQU())
	
	# Analytic solution in invariant form, matching integrate_ray's stepping convention.
	seg = S.z_interval(obj, ray)
	Δz = S.width(seg) / (ray.nz - 1)
	Δλ = (Δz / ν)
	Λtot = ray.nz * Δλ

	# u=0 => ν′=ν and invariants are simple.
	Jinv = j_mode / (ν^2)
	Ainv = a_mode * ν
	Iinv_mode(m) = (m[1] / m[2]) * (1 - exp(-m[2] * Λtot))
	Iinv_perp = Iinv_mode((Jinv.perp, Ainv.perp))
	Iinv_par = Iinv_mode((Jinv.par, Ainv.par))
	
	# Field Stokes basis uses +Q along e_perp, so Q_f = I_perp - I_par and U_f=0.
	S_f_inv = S.StokesIQU(Iinv_perp + Iinv_par, Iinv_perp - Iinv_par, 0.0)
	R_cf = S.stokes_QU_rotation(χ)
	R_fc = S.stokes_QU_rotation(-χ)
	QU_cam_inv = S.rotate_QU(R_fc, @swiz S_f_inv.QU)
	S_cam_inv = S.StokesIQU(S_f_inv.I, QU_cam_inv...)
	S_cam_expected = S_cam_inv * ν^3

	@test S_cam.I ≈ S_cam_expected.I rtol=1e-12
	@test S_cam.Q ≈ S_cam_expected.Q rtol=1e-12
	@test S_cam.U ≈ S_cam_expected.U rtol=1e-12
	@test abs(S_cam.U) > 0

	@testset "Depolarized limit matches scalar I" begin
		j0 = 0.8
		a0 = 1.3
		obj0 = ConstPolarizedSlab(z, u0, B′, S.ModePerpPar(j0/2, j0/2), S.ModePerpPar(a0, a0))
		I = S.render(ray, obj0, S.Intensity())
		S0 = S.render(ray, obj0, S.IntensityIQU())
		@test S0.I ≈ I
		@test S0.Q ≈ 0 atol=1e-12
		@test S0.U ≈ 0 atol=1e-12
	end
end


@testitem "Synchrotron slab polarization limits (thin/thick, velocity invariance)" begin
	import Synchray as S
	using Test

	p = 3.2
	model = S.IsotropicPowerLawElectrons(; p, Cj=0.7, Ca=0.3)
	ν = 2.0
	z = 0.0 .. 1.0
	ray = S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=S.photon_k(ν, SVector(0.0, 0.0, 1.0)), nz=80)

	Πj = S._Pi_j(p)
	Πα = S._Pi_a(p)
	Πthick = let
		Sperp = (1 + Πj) / (1 + Πα)
		Spar = (1 - Πj) / (1 - Πα)
		(Sperp - Spar) / (Sperp + Spar)
	end

	# Stokes convention (matches the standard astronomy/IAU usage in a fixed sky basis):
	# Q = I_x − I_y and U = I_{45°} − I_{135°}, where angles are measured in the camera
	# screen plane from +x toward +y (right-handed). Therefore for a polarization angle χ
	# relative to +x we expect (Q,U)/I = Π(cos(2χ), sin(2χ)).
	#
	# Here `χ` parameterizes the orientation of the field-aligned +Q axis (e_perp) relative
	# to the camera x-axis.
	# Choose B′ in the screen plane so the +Q axis (e_perp) aligns with camera x (U≈0).
	Bmag = 0.9
	B_align = SVector(0.0, Bmag, 0.0)
	
	# And a rotated field where the camera sees both Q and U.
	χ = 0.23
	B_rot = Bmag * SVector(-sin(χ), cos(χ), 0.0)

	u_rest = S.FourVelocity(SVector(0.0, 0.0, 0.0))
	u_move = S.FourVelocity(SVector(0.0, 0.0, 0.5))
	us = [u_rest, u_move]

	# Helper: scale ne0 so that the actual ray optical depth hits a target τ.
	with_tau(u0, B0, τ_target) = begin
		slab1 = S.UniformSynchrotronSlab(; z, u0, ne0=1.0, B0, model)
		τ1 = S.render(ray, slab1, S.OpticalDepth())
		@test τ1 > 0
		ne0 = τ_target / τ1
		S.UniformSynchrotronSlab(; z, u0, ne0, B0, model)
	end

	targets = (
		(τ=1e-4, Π=Πj),
		(τ=40.0, Π=Πthick),
	)

	@testset "Aligned field: U≈0, Q/I matches Π" begin
		for (; τ, Π) in targets
			vals = map(us) do u0
				slab = with_tau(u0, B_align, τ)
				S.render(ray, slab, S.IntensityIQU())
			end
			fracQ = map(s -> s.Q / s.I, vals)
			fracU = map(s -> s.U / s.I, vals)
			@test fracQ ≈ [Π, Π]  rtol=3e-3
			@test fracU ≈ [0.0, 0.0]  atol=1e-12
			@test fracQ[1] ≈ fracQ[2]
		end
	end

	@testset "Rotated field: direction matches 2χ law" begin
		for (; τ, Π) in targets
			vals = map(us) do u0
				slab = with_tau(u0, B_rot, τ)
				S.render(ray, slab, S.IntensityIQU())
			end
			fracQ = map(s -> s.Q / s.I, vals)
			fracU = map(s -> s.U / s.I, vals)
			Qexp = Π * cos(2χ)
			Uexp = Π * sin(2χ)
			@test fracQ ≈ [Qexp, Qexp] rtol=3e-3
			@test fracU ≈ [Uexp, Uexp] rtol=3e-3
			@test fracQ[1] ≈ fracQ[2]
			@test fracU[1] ≈ fracU[2]
		end
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

		# Choose B′ so that e_perp aligns with e1′ up to a 180° flip.
		(_, e_perp) = S.linear_polarization_basis_from_B(n′, SVector(0, 1, 0))
		R = S.stokes_QU_rotation(e1′, e2′, e_perp)
		@test R ≈ SMatrix{2,2}(I)
	end

	@testset "Known nonzero rotation in screen plane" begin
		ν = 2
		k = S.photon_k(ν, SVector(0, 0, 1))
		u = S.FourVelocity(SVector(0, 0, 0))
		(n′, e1′, e2′) = S.comoving_screen_basis(u, k)

		# Choose B′ so that e_perp (the +Q axis) is at a known angle χ from e1′.
		χ = 0.41
		B′ = sin(χ) * e1′ - cos(χ) * e2′
		(e_par, e_perp) = S.linear_polarization_basis_from_B(n′, B′)
		@test dot(e_par, B′) ≈ norm(B′)

		Rgeom = S.stokes_QU_rotation(e1′, e2′, e_perp)
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

		# Construct B′ with a nonzero component along n′, while e_perp aligns with e1′ (up to a 180° flip).
		B′ = -e2′ + 0.2 * n′
		(_, e_perp) = S.linear_polarization_basis_from_B(n′, B′)
		R = S.stokes_QU_rotation(e1′, e2′, e_perp)
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

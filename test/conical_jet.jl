@testitem "ConicalJet vs ConicalBKJet compatibility" begin
	import Synchray as S
	using Synchray: skip
	using RectiGrids

	axis = SVector(0.01, 0, 1)
	φj = 0.05
	s = 1e-3 .. 70
	s0 = 1
	ne0 = 5e6
	B0 = 1
	ne_exp = -2
	B_exp = -1

	xs = range(0 ± maximum(s) * tan(φj); length=10)
	cam = S.CameraZ(; xys=grid(SVector, xs, xs), nz=1024, ν=2, t=0)

	γs = (1.2, 10, 50)
	@testset for γ in γs
		speed_profile = η -> (S.gamma, γ)
		model = S.IsotropicPowerLawElectrons(p=2.5, Cj=1, Ca=1)

		bk = S.ConicalBKJet(; axis, φj, s, s0, ne0, B0, ne_exp, B_exp, speed_profile, model)
		cj = S.ConicalJet(; axis, φj, s, speed_profile,
			ne=S.PowerLawS(ne_exp; val0=ne0, s0=s0),
			B=S.BFieldSpec(S.PowerLawS(B_exp; val0=B0, s0=s0), S.ScalarField(), b -> S.FullyTangled(b)),
			model
		)
		cjp = S.prepare_for_computations(cj)

		@testset for what in (S.Intensity(), S.OpticalDepth(), S.SpectralIndex())
			img_bk = S.render(cam, bk, what)
			img_cj = S.render(cam, cj, what)
			img_cjp = S.render(cam, cjp, what)
			@test count(>(0), img_bk) > 40
			@test minimum(skip(isnan, img_bk)) < 0.9*maximum(skip(isnan, img_bk))
			@test img_cj ≈ img_bk  nans=true
			@test img_cjp ≈ img_bk  nans=true
		end
	end
end

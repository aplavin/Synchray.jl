@testitem "Unitful boundary API (ConicalJet)" begin
	import Synchray as S
	using Unitful, UnitfulAstro
	using RectiGrids

	jet = S.withunits(S.ConicalJet;
		axis=SVector(0, 0, 1),
		φj=2u"°",
		s=1e-3u"pc"..50u"pc",
		speed_profile=(η -> (S.beta, 0.0)),
		ne=S.PowerLawS(-2; val0=2u"cm^-3", s0=1u"pc"),
		B=S.BFieldSpec(S.PowerLawS(-1; val0=3u"Gauss", s0=1u"pc"), S.ScalarField(), b -> S.FullyTangled(b)),
		model=S.IsotropicPowerLawElectrons(; p=2.5),
	)

	cam = S.withunits(S.CameraZ;
		xys=grid(SVector, range(0.01u"pc"..0.1u"pc", 2), range(-0.001u"pc"..0.001u"pc", 2)),
		nz=20,
		ν=230u"GHz",
		t=0u"yr",
	)

	Iν_nou = S.render(cam, jet)
	@test all(>(0), Iν_nou)

	Iν = S.withunits(S.render, cam, jet, S.Intensity())
	@test unit(eltype(Iν)) == u"erg/s/cm^2/Hz/sr"
	@test ustrip.(Iν) == Iν_nou
	@test ustrip.(u"Jy/mas^2", Iν) ≈ [121 121; 0.028 0.028]  rtol=1e-2
end

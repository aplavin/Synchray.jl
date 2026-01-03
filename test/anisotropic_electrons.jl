@testitem "AnisotropicPowerLawElectrons: φ(θ) scaling" begin
	import Synchray as S
	using Test

	p = 2.5
	ν = 2
	n_e = 3
	b = S.SVector(0, 0, 5)

	θ = acos(0.5) # cosθ = 1/2
	n̂ = S.SVector(sin(θ), 0, cos(θ))
	k′ = S.FourFrequency(ν, (ν .* n̂)...)

	model_iso = S.IsotropicPowerLawElectrons(; p, Cj=1, Ca=1)
	model_an1 = S.AnisotropicPowerLawElectrons(; p, η=1, Cj=1, Ca=1)
	model_an = S.AnisotropicPowerLawElectrons(; p, η=0.5, Cj=1, Ca=1)

	(j_iso, α_iso) = S._synchrotron_coeffs(model_iso, n_e, b, k′)
	(j_an1, α_an1) = S._synchrotron_coeffs(model_an1, n_e, b, k′)
	(j_an, α_an) = S._synchrotron_coeffs(model_an, n_e, b, k′)

    @test j_iso > 0 && α_iso > 0

	@test j_an1 ≈ j_iso
	@test α_an1 ≈ α_iso

	cosθ = 0.5
	φ = (1 + (model_an.η - 1) * cosθ^2)^(-model_an.p / 2) / model_an.Pnorm
	@test j_an ≈ φ * j_iso
	@test α_an ≈ φ * α_iso
end

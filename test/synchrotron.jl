@testitem "Synchrotron slab scalings" begin
	import Synchray as S
	using Accessors

	L = 2
	ν = 2.0
	model = S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.0)
	ne0 = 1.3
	B0 = 0.9
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=ν, nz=4_000)

	cases = (
		S.FourVelocity(SVector(0, 0, 0)),
		S.FourVelocity(SVector(0, 0, 0.3)),
		S.FourVelocity(SVector(0.3, 0, 0.5)),
	)

	@testset for u0 in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		ν′ = ν / δ

		slab = S.UniformSynchrotronSlab(; z=0..L, u0, ne0, B0=S.FullyTangled(B0), model)

		@testset "optically thin" begin
			I_num = S.render(ray, slab)
			@test I_num > 0
			@test S.render(ray, slab, S.OpticalDepth()) ≈ 0

			k′ = S.photon_k(ν′, SVector(0, 0, 1))
			(j0, _) = S._synchrotron_coeffs(model, ne0, slab.B0, k′)
			I_exact = j0 * L * (δ^2)
			@test I_num ≈ I_exact rtol=2e-3

			@test S.render(ray, @set slab.ne0 *= 2) ≈ 2I_num rtol=2e-3
			@test S.render(ray, @set slab.B0 *= 2) ≈ I_num * (2^((model.p + 1) / 2)) rtol=2e-3
			@test S.render((@set S.photon_frequency(ray) *= 2), slab) ≈ I_num * (2^(-(model.p - 1) / 2)) rtol=2e-3

			αobs = S.render(ray, slab, S.SpectralIndex())
			@test αobs ≈ (-(model.p - 1) / 2) atol=2e-4
		end

		@testset "very optically thick" begin
			model_thick = @set model.Ca_ordered = 5e6
			slab_thick = @set slab.model = model_thick
			k′ = S.photon_k(ν′, SVector(0, 0, 1))
			(j_thick, α_thick) = S._synchrotron_coeffs(model_thick, ne0, slab.B0, k′)

			I_num_thick = S.render(ray, slab_thick)
			I_exact_thick = (j_thick / α_thick) * (δ^3)
			@test I_num_thick ≈ I_exact_thick rtol=2e-3
			@test S.render(ray, slab_thick, S.OpticalDepth()) ≈ (α_thick * (L / δ)) rtol=2e-3
			@test S.render(ray, slab_thick, S.SpectralIndex()) ≈ 5/2 atol=2e-4
			@test S.render((@set S.photon_frequency(ray) *= 2), slab_thick) ≈ I_num_thick * (2^(5/2)) rtol=2e-3

			@test S.render(ray, @set slab_thick.ne0 *= 2) ≈ I_num_thick rtol=1e-3
		end
	end
end


@testitem "Ordered-field synchrotron (minimal)" begin
	import Synchray as S
	using Accessors

	L = 2
	ν = 2.0
	p = 3.2
	model = S.IsotropicPowerLawElectrons(; p, Cj=0.7, Ca=0)
	ne0 = 1.3
	B0 = 0.9

	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=ν, nz=4_000)
	u0 = S.FourVelocity(SVector(0, 0, 0))
	k′ = S.photon_k(ν, SVector(0, 0, 1))

	@testset "B parallel to photon => zero emissivity" begin
		slab = S.UniformSynchrotronSlab(; z=0..L, u0, ne0, B0=SVector(0, 0, B0), model)
		I = S.render(ray, slab)
		@test I ≈ 0
		(j, α) = S._synchrotron_coeffs(model, ne0, slab.B0, k′)
		@test j ≈ 0
		@test α ≈ 0
	end

	@testset "B ⟂ photon => thin-slab analytic limit" begin
		slab = S.UniformSynchrotronSlab(; z=0..L, u0, ne0, B0=SVector(B0, 0, 0), model)
		I_num = S.render(ray, slab)
		(j, _) = S._synchrotron_coeffs(model, ne0, slab.B0, k′)
		I_exact = j * L
		@test I_num ≈ I_exact  rtol=1e-3
	end

	@testset "sin(theta) scaling" begin
		slab_perp = S.UniformSynchrotronSlab(; z=0..L, u0, ne0, B0=SVector(B0, 0, 0), model)
		I_perp = S.render(ray, slab_perp)

		slab_45 = @set slab_perp.B0 = B0 * normalize(SVector(1, 0, 1))

		expected = (1 / sqrt(2))^((p + 1) / 2)
		@test (S.render(ray, slab_45) / I_perp) ≈ expected
	end
end


@testitem "TangledOrderedMixture magnetic field" begin
	import Synchray as S

	p = 3.2
	ne0 = 1.3
	B0 = 0.9
	ν = 2.0
	k′ = S.photon_k(ν, SVector(0.0, 0.0, 1.0))

	model = S.IsotropicPowerLawElectrons(; p)

	# A generic non-perpendicular field direction.
	b = B0 .* normalize(SVector(1.0, 0.0, 2.0))
	@testset "limits" begin
		(j_t, α_t) = S._synchrotron_coeffs(model, ne0, S.FullyTangled(B0), k′)
		(j_0, α_0) = S._synchrotron_coeffs(model, ne0, S.TangledOrderedMixture(b; kappa=0.0), k′)
		(j_inf, α_inf) = S._synchrotron_coeffs(model, ne0, S.TangledOrderedMixture(b; kappa=Inf), k′)
		(j_ord, α_ord) = S._synchrotron_coeffs(model, ne0, b, k′)

		@test j_0 ≈ j_t
		@test α_0 ≈ α_t
		@test j_inf ≈ j_ord
		@test α_inf ≈ α_ord
	end

	@testset "intermediate kappa matches interpolation" begin
		κ = 2.0
		f = κ / (1 + κ)
		(j_mid, α_mid) = S._synchrotron_coeffs(model, ne0, S.TangledOrderedMixture(b; kappa=κ), k′)

		B = norm(b)
		(j_perp, α_perp) = S._synchrotron_coeffs(model, ne0, SVector(B, 0.0, 0.0), k′)

		μ = b.z / B
		sinθ2 = clamp(1 - μ^2, 0, 1)
		qj = (p + 1) / 2
		qa = (p + 2) / 2
		sinpow_j = sinθ2^(qj / 2)
		sinpow_a = sinθ2^(qa / 2)
		Aj = (1 - f) * model.sinavg_j + f * sinpow_j
		Aa = (1 - f) * model.sinavg_a + f * sinpow_a

		@test j_mid ≈ j_perp * Aj
		@test α_mid ≈ α_perp * Aa

		(j_0, α_0) = S._synchrotron_coeffs(model, ne0, S.TangledOrderedMixture(b; kappa=0.0), k′)
		(j_inf, α_inf) = S._synchrotron_coeffs(model, ne0, S.TangledOrderedMixture(b; kappa=Inf), k′)
		@test min(j_0, j_inf) ≤ j_mid ≤ max(j_0, j_inf)
		@test min(α_0, α_inf) ≤ α_mid ≤ max(α_0, α_inf)
	end
end


@testitem "Ordered vs tangled consistency" begin
	import Synchray as S

	avg_sin_pow(q) = sqrt(pi) * S.SpecialFunctions.gamma((q + 2) / 2) / (2 * S.SpecialFunctions.gamma((q + 3) / 2))

	p = 3.2
	ne0 = 1.3
	B0 = 0.9
	ν = 2.0
	k′ = S.photon_k(ν, SVector(0.0, 0.0, 1.0))

	model = S.IsotropicPowerLawElectrons(; p)

	(j_t, α_t) = S._synchrotron_coeffs(model, ne0, S.FullyTangled(B0), k′)
	(j_o_perp, α_o_perp) = S._synchrotron_coeffs(model, ne0, SVector(B0, 0.0, 0.0), k′)

	@testset "analytic ratio vs orthogonal field" begin
		qj = (p + 1) / 2
		qa = (p + 2) / 2
		@test (j_t / j_o_perp) ≈ avg_sin_pow(qj) rtol=2e-12
		@test (α_t / α_o_perp) ≈ avg_sin_pow(qa) rtol=2e-12
	end

	@testset "numerical angle average matches tangled" begin
		θs = range(0.0, pi; length=100)
		ws = sin.(θs)

		js = map(θ -> begin
			b = B0 .* SVector(sin(θ), 0.0, cos(θ))
			first(S._synchrotron_coeffs(model, ne0, b, k′))
		end, θs)
		αs = map(θ -> begin
			b = B0 .* SVector(sin(θ), 0.0, cos(θ))
			last(S._synchrotron_coeffs(model, ne0, b, k′))
		end, θs)

		j_avg = sum(js .* ws) / sum(ws)
		α_avg = sum(αs .* ws) / sum(ws)

		@test j_avg ≈ j_t rtol=2e-3
		@test α_avg ≈ α_t rtol=2e-3
	end
end

@testitem "AnisotropicPowerLawElectrons: φ(θ) scaling" begin
	import Synchray as S
	using Test

	p = 2.5
	ν = 2
	n_e = 3
	b = SVector(0, 0, 5)

	θ = acos(0.5) # cosθ = 1/2
	n̂ = SVector(sin(θ), 0, cos(θ))
	k′ = S.FourFrequency(ν, (ν .* n̂)...)

	model_iso = S.IsotropicPowerLawElectrons(; p, Cj=1, Ca=1)
	model_an1 = S.AnisotropicPowerLawElectrons(; p, η=1, Cj=1, Ca=1)
	model_an = S.AnisotropicPowerLawElectrons(; p, η=0.5, Cj=1, Ca=1)

	(j_iso, α_iso) = S._synchrotron_coeffs(model_iso, n_e, b, k′)
	@test j_iso > 0 && α_iso > 0

	(j_an1, α_an1) = S._synchrotron_coeffs(model_an1, n_e, b, k′)
	@test j_an1 ≈ j_iso
	@test α_an1 ≈ α_iso

	(j_an, α_an) = S._synchrotron_coeffs(model_an, n_e, b, k′)
	cosθ = 0.5
	φ = (1 + (model_an.η - 1) * cosθ^2)^(-model_an.p / 2) / model_an.Pnorm
	@test j_an ≈ φ * j_iso
	@test α_an ≈ φ * α_iso
end

@testitem "prepare_for_computations preserves synchrotron coeffs" begin
	import Synchray as S
	using Test

	ne0 = 1.3
	ν = 2
	k′ = S.photon_k(ν, SVector(0, 0, 1))

	models = (
		S.IsotropicPowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.3),
		S.AnisotropicPowerLawElectrons(; p=2.5, η=0.5, Cj=1, Ca=1),
	)

	b = SVector(1, 0, 2)
	B_ordered = (SVector(3, 0, 0), SVector(0, 0, 3))
	B_tangled = (
		S.FullyTangled(3),
		S.TangledOrderedMixture(b; kappa=0),
		S.TangledOrderedMixture(b; kappa=2),
		S.TangledOrderedMixture(b; kappa=Inf),
	)

	@testset for model in models
		prepared = S.prepare_for_computations(model)
		Bs = model isa S.AnisotropicPowerLawElectrons ? B_ordered : (B_ordered..., B_tangled...)

		@testset for B in Bs
			(j1, α1) = S._synchrotron_coeffs(model, ne0, B, k′)
			(j2, α2) = S._synchrotron_coeffs(prepared, ne0, B, k′)

			B == SVector(0,0,3) || @test j1 > 0 && α1 > 0
			@test j1 ≈ j2 rtol=0.5e-3
			@test α1 ≈ α2 rtol=0.5e-3
		end
	end
end

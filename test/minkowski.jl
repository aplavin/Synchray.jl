@testitem "Minkowski + Doppler conventions" begin
	import Synchray as S

	@testset "arithmetics" begin
		a = S.FourPosition(1, 2, 3, 4)
		b = S.FourPosition(0.5, -1.0, 2.0, 3.0)
		@test a + b === S.FourPosition(1.5, 1.0, 5.0, 7.0)
		@test a - b === S.FourPosition(0.5, 3.0, 1.0, 1.0)
		@test 2a === S.FourPosition(2, 4, 6, 8)
	end

	@testset "minkowski_dot basics" begin
		a = S.FourPosition(2.0, 1.0, 3.0, -4.0)
		b = S.FourPosition(-1.0, 2.0, 0.5, 7.0)
		@test S.minkowski_dot(a, b) ≈ S.minkowski_dot(b, a)
		@test S.minkowski_dot(a, a) ≈ -(a.t^2) + a.x^2 + a.y^2 + a.z^2
	end

	@testset "FourVelocity normalization" begin
		β = SVector(0.3, -0.4, 0.1)
		u = S.FourVelocity(β)
		@test S.minkowski_dot(u, u) ≈ -1
		@test u.t ≈ inv(sqrt(1 - dot(β, β)))
		@test S.beta(u) ≈ β
		@test S.gamma(u) ≈ u.t
		@test S.gamma(β) ≈ u.t
	end

	@testset "photon_k is null" begin
		ν = 2.5
		n = normalize(SVector(0.2, -0.3, 0.7))
		k = S.photon_k(ν, n)
		@test S.minkowski_dot(k, k) ≈ 0  atol=√eps(Float64)
	end

	@testset "Doppler boosting δ convention" begin
		# Convention locked in by this test:
		# doppler_factor(u, n) == δ = ν_obs / ν_comoving
		# For fast motion toward the observer along the photon direction (β ⋅ n > 0): δ >> 1.
		β = 0.99
		u_toward = S.FourVelocity(SVector(0.0, 0.0, β))
		n = SVector(0.0, 0.0, 1.0)

		δ_toward = S.doppler_factor(u_toward, n)
		@test δ_toward > 1
		@test δ_toward > 10

		u_away = S.FourVelocity(SVector(0.0, 0.0, -β))
		δ_away = S.doppler_factor(u_away, n)
		@test δ_away < 1
	end

	@testset "Lorentz boost (rest → lab)" begin
		β = 0.6
		u = S.FourVelocity(SVector(0, 0, β))
		γ = u.t

		Λ = S.lorentz_boost_matrix(u)
		v = S.FourPosition(1, 2, 0.5, -1)
		v′ = S.lorentz_boost(u, v)::S.FourPosition
		@test v′.x ≈ v.x
		@test v′.y ≈ v.y
		@test v′.t ≈ γ * (v.t + β * v.z)
		@test v′.z ≈ γ * (v.z + β * v.t)
		@test S.minkowski_dot(v′, v′) ≈ S.minkowski_dot(v, v)
		@test Λ * v ≈ v′

		u0 = S.FourVelocity(1, 0, 0, 0)
		@test S.lorentz_boost(u, u0)::S.FourVelocity ≈ u
		@test Λ * u0 ≈ u

		uid = S.FourVelocity(SVector(0, 0, 0))
		Λid = S.lorentz_boost_matrix(uid)
		@test Λid * v ≈ v
		@test S.lorentz_boost(uid, v) ≈ v

		u_neg = S.FourVelocity(u.t, -@swiz u.xyz)
		@test S.lorentz_boost_matrix(u_neg) ≈ inv(Λ)
		@test S.lorentz_boost(u_neg, v′) ≈ v
	end

	@testset "lorentz_unboost is inverse" begin
		β = SVector(0.3, -0.4, 0.1)
		u = S.FourVelocity(β)

		vpos = S.FourPosition(0.7, 1.2, -0.5, 3.4)
		@test S.lorentz_unboost(u, S.lorentz_boost(u, vpos)) ≈ vpos
		@test S.lorentz_boost(u, S.lorentz_unboost(u, vpos)) ≈ vpos

		n = normalize(SVector(0.2, 0.1, 0.7))
		k = S.photon_k(2.3, n)
		@test S.lorentz_unboost(u, S.lorentz_boost(u, k)) ≈ k
		@test S.lorentz_boost(u, S.lorentz_unboost(u, k)) ≈ k

		uid = S.FourVelocity(SVector(0, 0, 0))
		@test S.lorentz_unboost(uid, vpos) === vpos
		@test S.lorentz_unboost(uid, k) === k
	end
end

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
		@test u.t ≈ inv(√(1 - dot(β, β)))
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

@testitem "proper_velocity (celerity)" begin
    import Synchray as S

    β = S.SVector(0.1, -0.2, 0.3)
    u = S.FourVelocity(β)                                   # existing ctor (γ, γβ⃗)
    w = S.proper_velocity(u)
    @test w ≈ u.t * β
    @test S.beta(u) ≈ w / u.t
    # cancellation-free reconstruction round-trips:
    @test S.proper_velocity(S.construct(S.FourVelocity, S.proper_velocity => w)) ≈ w
    @test S.construct(S.FourVelocity, S.proper_velocity => w) ≈ S.FourVelocity(β)
    # high-γ: build from celerity stays well-conditioned where β→1 is ill-represented
    γ = 50.0; βmag = sqrt(1 - 1/γ^2); whi = (γ*βmag) * S.SVector(1.0, 0.0, 0.0)
    uhi = S.construct(S.FourVelocity, S.proper_velocity => whi)
    @test uhi.t ≈ γ
end

@testitem "Float32 lorentz_unboost null-vector stability" begin
    import Synchray as S

    # For a null FourFrequency, |k.xyz| = k.t; unboost must preserve this exactly.
    # In Float32 the generic boost suffers catastrophic cancellation when γ is large
    # and the photon is roughly aligned with β. Assert Float32 and Float64 agree to
    # ~Float32 precision.

    @testset "γ=$γ, Δθ=$(Δθ)°" for γ in (5.0, 30.0, 46.0, 100.0),
                                    Δθ in (0.0, 1.0, 2.0, 6.0, 90.0)
        β = sqrt(1 - 1/γ^2)
        axis = SVector(sind(Δθ), 0.0, cosd(Δθ))

        u64 = S.FourVelocity(γ, γ*β*axis...)
        k64 = S.FourFrequency(1.0, 0.0, 0.0, 1.0)            # photon along +z

        u32 = S.to_float_type(Float32, u64)
        k32 = S.to_float_type(Float32, k64)

        kp64 = S.lorentz_unboost(u64, k64)
        kp32 = S.lorentz_unboost(u32, k32)

        # Null condition preserved → |xyz'|² == t'². rtol tighter than the default `≈`
        # (√eps(Float32) ≈ 3.5e-4) to catch the γ=46 cancellation regime.
        x64 = SVector(kp64.x, kp64.y, kp64.z)
        x32 = SVector(kp32.x, kp32.y, kp32.z)
        @test dot(x64, x64) ≈ kp64.t^2  rtol=1e-12
        @test dot(x32, x32) ≈ kp32.t^2  rtol=1e-4

        # Float32 result agrees with Float64 reference to ~Float32 precision.
        @test kp32.t ≈ Float32(kp64.t)            rtol=1e-4
        @test x32    ≈ SVector{3,Float32}(x64)    rtol=1e-4
    end
end

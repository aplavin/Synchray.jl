@testitem "CombinedMedium" begin
    import Synchray as S
    using RectiGrids

    # Three non-overlapping stationary ellipsoids at increasing z (back-to-front)
    ell_A = S.MovingUniformEllipsoid(
        center=S.FourPosition(0, 0, 0, 3.0),
        sizes=SVector(1.0, 1.0, 0.5),
        u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)),
        jν=0.7, αν=1.3,
    )
    ell_B = S.MovingUniformEllipsoid(
        center=S.FourPosition(0, 0, 0, 6.0),
        sizes=SVector(0.8, 0.8, 0.6),
        u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)),
        jν=1.2, αν=0.5,
    )
    ell_C = S.MovingUniformEllipsoid(
        center=S.FourPosition(0, 0, 0, 10.0),
        sizes=SVector(0.7, 0.7, 0.4),
        u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)),
        jν=0.9, αν=0.8,
    )

    xs = range(-1.5..1.5, 32)
    cam = S.CameraZ(; xys=grid(SVector, xs, xs), nz=512, ν=2, t=0)

    @testset "single-object passthrough" begin
        cm1 = S.CombinedMedium(ell_A)
        @test S.render(cam, cm1) ≈ S.render(cam, ell_A)
        @test S.render(cam, cm1, S.OpticalDepth()) ≈ S.render(cam, ell_A, S.OpticalDepth())
    end

    @testset "two objects — all modes" begin
        cm = S.CombinedMedium(ell_A, ell_B)  # A is back (lower z), B is front

        # OpticalDepth: additive for non-overlapping objects
        img_τ_A = S.render(cam, ell_A, S.OpticalDepth())
        img_τ_B = S.render(cam, ell_B, S.OpticalDepth())
        img_τ_cm = S.render(cam, cm, S.OpticalDepth())
        @test img_τ_cm ≈ img_τ_A .+ img_τ_B rtol=2e-3

        # Intensity: back object A attenuated by front object B
        img_I_A = S.render(cam, ell_A)
        img_I_B = S.render(cam, ell_B)
        img_I_cm = S.render(cam, cm)
        @test img_I_cm ≈ img_I_A .* exp.(.-img_τ_B) .+ img_I_B rtol=2e-3

        # SpectralIndex: frequency-independent j, α → α_spec = 0
        img_α = S.render(cam, cm, S.SpectralIndex())
        @test all(x -> isnan(x) || abs(x) < 1e-5, img_α)

        # Tuple combos
        img_I_τ = S.render(cam, cm, (S.Intensity(), S.OpticalDepth()))
        @test getindex.(img_I_τ, 1) ≈ img_I_cm
        @test getindex.(img_I_τ, 2) ≈ img_τ_cm

        img_I_α = S.render(cam, cm, (S.Intensity(), S.SpectralIndex()))
        @test getindex.(img_I_α, 1) ≈ img_I_cm
    end

    @testset "three objects" begin
        cm = S.CombinedMedium(ell_A, ell_B, ell_C)  # back-to-front order

        # τ is additive
        img_τ_A = S.render(cam, ell_A, S.OpticalDepth())
        img_τ_B = S.render(cam, ell_B, S.OpticalDepth())
        img_τ_C = S.render(cam, ell_C, S.OpticalDepth())
        img_τ_cm = S.render(cam, cm, S.OpticalDepth())
        @test img_τ_cm ≈ img_τ_A .+ img_τ_B .+ img_τ_C rtol=2e-3

        # Intensity: layered attenuation A→B→C
        img_I_A = S.render(cam, ell_A)
        img_I_B = S.render(cam, ell_B)
        img_I_C = S.render(cam, ell_C)
        img_I_cm = S.render(cam, cm)
        @test img_I_cm ≈ img_I_A .* exp.(.-img_τ_B .- img_τ_C) .+ img_I_B .* exp.(.-img_τ_C) .+ img_I_C rtol=2e-3

        @test isbits(cm)
    end

    @testset "missed rays" begin
        xs_far = range(10.0..12.0, 8)
        cam_far = S.CameraZ(; xys=grid(SVector, xs_far, xs_far), nz=128, ν=2, t=0)
        cm = S.CombinedMedium(ell_A, ell_B)
        @test all(==(0), S.render(cam_far, cm))
        @test all(==(0), S.render(cam_far, cm, S.OpticalDepth()))
        @test all(isnan, S.render(cam_far, cm, S.SpectralIndex()))
    end

    @testset "ray hitting only one object" begin
        # ell_A has radius 1.0 in x,y; ell_B has 0.8 — camera at x≈0.9 hits A but misses B
        xs_partial = range(0.85..0.95, 4)
        ys_partial = range(-0.05..0.05, 4)
        cam_partial = S.CameraZ(; xys=grid(SVector, xs_partial, ys_partial), nz=512, ν=2, t=0)
        cm = S.CombinedMedium(ell_A, ell_B)
        @test S.render(cam_partial, cm) ≈ S.render(cam_partial, ell_A) rtol=1e-4
        @test S.render(cam_partial, cm, S.OpticalDepth()) ≈ S.render(cam_partial, ell_A, S.OpticalDepth()) rtol=1e-4
    end

    @testset "moving ellipsoids" begin
        ell_mov1 = S.MovingUniformEllipsoid(
            center=S.FourPosition(0, 0, 0, 2.0),
            sizes=SVector(0.5, 0.5, 0.5),
            u0=S.FourVelocity(SVector(0.0, 0.0, 0.3)),
            jν=1.0, αν=0.5,
        )
        ell_mov2 = S.MovingUniformEllipsoid(
            center=S.FourPosition(0, 0, 0, 5.0),
            sizes=SVector(0.5, 0.5, 0.5),
            u0=S.FourVelocity(SVector(0.0, 0.0, -0.2)),
            jν=0.8, αν=0.3,
        )
        cm = S.CombinedMedium(ell_mov1, ell_mov2)
        @test sum(S.render(cam, cm)) > 0
    end

    @testset "isbits" begin
        @test isbits(S.CombinedMedium(ell_A, ell_B))
        @test isbits(S.CombinedMedium(ell_A, ell_B, ell_C))
    end

    @testset "vector of objects matches tuple" begin
        cm_tuple = S.CombinedMedium(ell_A, ell_B, ell_C)
        cm_vec = S.CombinedMedium([ell_A, ell_B, ell_C])

        @test S.render(cam, cm_vec) ≈ S.render(cam, cm_tuple)
        @test S.render(cam, cm_vec, S.OpticalDepth()) ≈ S.render(cam, cm_tuple, S.OpticalDepth())
    end
end

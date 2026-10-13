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
        @test sum(S.render(cam, cm)) ≈ 117.6 rtol=1e-2
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

@testitem "CombinedMedium: pointwise medium and per-component steps" begin
    import Synchray as S

    fixed(geometry, j, α) = S.EmissionRegion(; geometry, emission=S.FixedEmission(S=j / α, α=α))
    ball(R, j, α) = fixed(S.Geometries.Ball(S.Geometries.StationaryAtOrigin(), R), j, α)
    # different component types: stationary ball, and a ball on a moving worldline around z = 8
    back = ball(5.0, 0.7, 1.3)
    front = fixed(S.Geometries.Ball(; center=S.Geometries.InertialWorldline(x0=S.FourPosition(0.0, 0.0, 0.0, 8.0), u=S.FourVelocity(SVector(0.0, 0.0, 0.3))), size=1.0), 1.2, 0.5)
    inner = ball(3.0, 0.2, 0.4)
    cm = S.CombinedMedium(back, front)
    overlapping = (S.CombinedMedium(back, inner), S.CombinedMedium([back, inner]))

    x(z) = S.FourPosition(0.0, 0.0, 0.0, z)
    k′ = S.FourFrequency(2.0, 0.0, 0.0, 2.0)
    @testset "pointwise: first containing component" begin
        @test S.emissivity_absorption(cm, x(2.0), k′) == (0.7, 1.3)
        @test S.emissivity_absorption(cm, x(8.0), k′) == (1.2, 0.5)
        @test S.four_velocity(cm, x(8.0)) == S.four_velocity(front, x(8.0))
        # outside all components: vacuum, in the first component's frame
        @test S.emissivity_absorption(cm, x(6.5), k′) == (0.0, 0.0)
        @test S.four_velocity(cm, x(6.5)) == S.four_velocity(back, x(6.5))
        @testset for cmo in overlapping
            @test S.emissivity_absorption(cmo, x(2.0), k′) == (0.7, 1.3)
        end
    end

    st = S.RayState(S.RayZ(; x0=S.FourPosition(0.0, 0.0, 0.0, 0.0), k=2.0, nz=1))
    Δλ = 0.05
    step(acc, obj, z) = S._integrate_ray_step(acc, obj, x(z), st, Δλ)
    @testset "step: each containing component" begin
        acc = S.AccValue{S.Intensity}(0.3)
        # Δτ = 1.3·2·0.05: 0.3·e^(−Δτ) + (0.7/1.3)(1 − e^(−Δτ))
        @test step(acc, cm, 2.0).value ≈ 0.3290695510881738
        @testset for what in (S.Intensity, S.OpticalDepth, Tuple{S.Intensity, S.OpticalDepth}, Tuple{S.Intensity, S.SpectralIndex})
            acc = S._init_acc(what, 2.0)
            @test step(acc, cm, 2.0) === step(acc, back, 2.0)
            @test step(acc, cm, 8.0) === step(acc, front, 8.0)
            @test step(acc, cm, 6.5) === acc
            @testset for cmo in overlapping
                @test step(acc, cmo, 2.0) === step(step(acc, back, 2.0), inner, 2.0)
            end
        end
    end

    @testset "mixed components: type-stable, allocation-free" begin
        sync(geometry, Bdir) = S.EmissionRegion(; geometry,
            emission=S.SynchrotronEmission(
                ne=S.Profiles.Radial(r -> exp(-(r / 4)^2)),
                B=S.BFieldSpec(S.Profiles.Constant(1.0), S.Directions.Scalar(), Bdir),
                electrons=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1e-2))) |> S.prepare_for_computations
        cm_sync = S.CombinedMedium(
            sync(S.Geometries.Ball(S.Geometries.StationaryAtOrigin(), 5.0), b -> SVector(zero(b), zero(b), b)),
            sync(S.Geometries.Ball(; center=S.Geometries.InertialWorldline(x0=S.FourPosition(0.0, 0.0, 0.0, 8.0), u=S.FourVelocity(SVector(0.3, 0.0, 0.1))), size=1.0), b -> S.FullyTangled(b)))
        allocs(acc, obj, x4, st, Δλ) = (S._integrate_ray_step(acc, obj, x4, st, Δλ); @allocated S._integrate_ray_step(acc, obj, x4, st, Δλ))
        @testset for what in (S.Intensity, S.OpticalDepth, Tuple{S.Intensity, S.SpectralIndex}), z in (2.0, 6.5, 8.0)
            acc = S._init_acc(what, 2.0)
            @test (@inferred S._integrate_ray_step(acc, cm_sync, x(z), st, Δλ)) isa typeof(acc)
            @test allocs(acc, cm_sync, x(z), st, Δλ) == 0
        end
        # polarized: ordered fields of different types
        cm_pol = S.CombinedMedium(cm_sync.objects[1], sync(S.Geometries.Ball(; center=S.Geometries.InertialWorldline(x0=S.FourPosition(0.0, 0.0, 0.0, 8.0), u=S.FourVelocity(SVector(0.3, 0.0, 0.1))), size=1.0), b -> SVector(zero(b), b / 2, b)))
        @testset for z in (2.0, 6.5, 8.0)
            acc = S._init_acc(S.IntensityIQU, 2.0)
            @test (@inferred S._integrate_ray_step(acc, cm_pol, x(z), st, Δλ)) isa typeof(acc)
            @test allocs(acc, cm_pol, x(z), st, Δλ) == 0
        end
    end
end

@testitem "CombinedMedium: rays missing the first object" begin
    import Synchray as S
    using RectiGrids

    # Counterjet (back) + jet (front); the camera only sees the jet, so every ray misses the first object.
    jet(dir, electrons, field_dir, wrap) = S.EmissionRegion(;
        geometry = S.Geometries.Conical(; axis = dir * SVector(sin(0.1), 0, cos(0.1)), φj = 0.05, z = 0.1..50.0),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Constant(3.0)),
        emission = S.SynchrotronEmission(;
            ne = S.Profiles.Axial(S.PowerLaw(-2.0; val0 = 1.0, s0 = 1.0)),
            B = S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1.0; val0 = 1.0, s0 = 1.0)), field_dir, wrap),
            electrons),
    ) |> S.prepare_for_computations
    cam = S.CameraZ(; xys = grid(SVector, range(0.5..4, 12), range(-0.2..0.2, 5)), nz = 64, ν = 1.0, t = 0.0)

    iso = S.IsotropicPowerLawElectrons(; p = 2.5, Cj = 1.0, Ca = 1.0)
    aniso = S.AnisotropicPowerLawElectrons(; p = 2.5, η = 2.0, Cj = 1.0, Ca = 1.0)
    helical = S.Directions.HelicalAT(0.5)
    configs = [
        (iso,   S.Directions.Scalar(), b -> S.FullyTangled(b)),
        (iso,   helical, b -> S.TangledOrderedMixture(b; kappa = 1.0)),
        (iso,   helical, identity),
        (aniso, helical, b -> S.TangledOrderedMixture(b; kappa = 1.0)),
        (aniso, helical, identity),
    ]
    @testset for (electrons, field_dir, wrap) in configs, what in (S.IntensityIQU(), S.SpectralIndex(), S.OpticalDepth())
        front = jet(+1, electrons, field_dir, wrap)
        alone = S.render(cam, front, what)
        @test count(x -> all(isfinite, x), alone) > 10
        @test isequal(S.render(cam, S.CombinedMedium(jet(-1, electrons, field_dir, wrap), front), what), alone)
    end

    # Root cause: at zero field (e.g. outside a region) coefficients must be exactly 0 with zero ν-derivative.
    k′ = S.photon_k(S.ForwardDiff.Dual(1.0, 1.0), SVector(0.0, 0.0, 1.0))
    b0 = SVector(0.0, 0.0, 0.0)
    fields = [(iso, S.FullyTangled(0.0)), (iso, b0), (iso, S.TangledOrderedMixture(b0; kappa = 1.0)),
              (aniso, b0), (aniso, S.TangledOrderedMixture(b0; kappa = 1.0))]
    @testset for (electrons, field) in fields
        (j, α) = S._synchrotron_coeffs(S.prepare_for_computations(electrons), 1.0, field, k′)
        @test (S.ForwardDiff.value(j), S.ForwardDiff.value(α)) == (0, 0)
        @test (S.ForwardDiff.partials(j)[1], S.ForwardDiff.partials(α)[1]) == (0, 0)
    end
end

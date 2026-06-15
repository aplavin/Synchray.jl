include("helpers.jl")


function bk_jet_1_knot_snapshots_image(; render_field, log_path, suffix="")
    log_section!(log_path, "bk_jet_1_knot_snapshots_image")

    φj = 4u"°"
    θ = 3 * φj  # "small viewing angle", same as in bk_jet_image
    axis = SVector(sin(θ), 0, cos(θ))

    knot = let
        x_c0 = S.FourPosition(0, 0.1 * axis)
        u = S.FourVelocity(0.995 * axis)
        sizing = S.Patterns.CrossSectionSizing(0.3, 0.5)
        S.Patterns.EllipsoidalKnot(; x_c0, u, sizing, profile=S.Patterns.GaussianBump(100))
    end

    jet = S.EmissionRegion(
        geometry = S.Geometries.Conical(; axis, φj, z = 1e-3 .. 50.0),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.beta, S.Profiles.Constant(0.995)),
        emission = S.SynchrotronEmission(
            ne = S.Profiles.Modified(S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)), knot),
            B = S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
            electrons = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=1),
        ),
    ) |> S.prepare_for_computations

    ts = [0, 0.1, 0.25]
    whats = [
        (name="Intensity", what=S.Intensity(), kwargs=p -> (; colormap=:magma, colorscale=SymLog(1e-4 * p))),
        (name="α", what=S.SpectralIndex(), kwargs=_ -> (; colormap=:balance, colorrange=(-3, 3))),
    ]

    fig = Figure()
    for (r, w) in enumerate(whats), (c, t) in enumerate(ts)
        pos = fig[r, c]
        Axis(pos[1,1]; title="$(w.name), t=$(t)", aspect=DataAspect(), width=300, height=300)
        img = render_field(jet; extent=3, ν=1, t, what=w.what)
        p = w.what isa S.Intensity ? maximum(img) : 1
        plt = heatmap!(img; w.kwargs(p)...)
        Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_1_knot_snapshots$(suffix).png"), fig)
    fig
end

function precessing_nozzle_snapshots_image(; render_field, log_path, suffix="")
    log_section!(log_path, "precessing_nozzle_snapshots_image")

    φj = 4u"°"
    θ = 3 * φj
    axis = SVector(sin(θ), 0, cos(θ))

    nozzle = S.Patterns.PrecessingNozzle(
        θ_precession = 0.02,
        θ_nozzle = 0.005,
        period = 50.0,
        β_flow = 0.995,
        profile = S.Patterns.TophatBump(100.0),
    )

    jet = S.EmissionRegion(
        geometry = S.Geometries.Conical(; axis, φj, z = 1e-3 .. 50.0),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.beta, S.Profiles.Constant(0.995)),
        emission = S.SynchrotronEmission(;
            ne = S.Profiles.Modified(S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)), nozzle),
            B = S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
            electrons = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=1),
        ),
    ) |> S.prepare_for_computations

    ts = [0, 12.5, 25.0]  # period=50, so 0, T/4, T/2
    whats = [
        (name="Intensity", what=S.Intensity(), kwargs=p -> (; colormap=:magma, colorscale=SymLog(1e-4 * p))),
        (name="α", what=S.SpectralIndex(), kwargs=_ -> (; colormap=:balance, colorrange=(-3, 3))),
    ]

    fig = Figure()
    for (r, w) in enumerate(whats), (c, t) in enumerate(ts)
        pos = fig[r, c]
        Axis(pos[1,1]; title="$(w.name), t=$(t)", aspect=DataAspect(), width=300, height=300)
        img = render_field(jet; extent=3, ν=1, t, what=w.what)
        p = w.what isa S.Intensity ? maximum(img) : 1
        plt = heatmap!(img; w.kwargs(p)...)
        Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "precessing_nozzle_snapshots$(suffix).png"), fig)
    fig
end


function main(; render_field, log_path, suffix="")
    bk_jet_1_knot_snapshots_image(; render_field, log_path, suffix)
    precessing_nozzle_snapshots_image(; render_field, log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
    logs = init_script(@__FILE__)
    main(; render_field=make_render_cpu(logs.log_cpu), log_path=logs.log_cpu)
    if HAS_METAL
        main(; render_field=make_render_metal(logs.log_metal), log_path=logs.log_metal, suffix="_metal")
    end
end

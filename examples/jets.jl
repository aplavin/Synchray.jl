include("helpers.jl")


function bk_jet_image(; render_field, log_path, suffix="", is_cpu=true)
    log_section!(log_path, "bk_jet_image")

    φj = 4u"°"
    axis_for_viewing_angle(θ) = SVector(sin(θ), 0.0, cos(θ))

    jet0 = S.EmissionRegion(
        geometry = Geometries.Conical(; axis=SVector(1, 0, 0), φj, z=1e-3 .. 50),
        velocity = S.VelocitySpec(Directions.Radial(), S.gamma, Profiles.Transverse(Profiles.LinearInterp(((0.3, 10.0), (0.8, 8.0))))),
        emission = S.SynchrotronEmission(
            ne = Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
            B = S.BFieldSpec(Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), Directions.Scalar(), b->S.FullyTangled(b)),
            electrons = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1.0, Ca=0.1),
        ),
    ) |> S.prepare_for_computations

    views = (
        (name="inside cone", θ=0.5 * φj),
        (name="small angle", θ=3 * φj),
        (name="large angle", θ=45u"°"),
        (name="counter-jet-like", θ=180u"°" - 3 * φj),
    )
    whats = [
        (what=S.Intensity(), kwargs=p -> (;colormap=:magma, colorscale=SymLog(1e-4*p))),
        (what=S.SpectralIndex(), kwargs=p -> (;colormap=:balance, colorrange=(-3, 3))),
        (what=S.OpticalDepth(), kwargs=p -> (;colormap=:viridis, colorrange=(0, min(6, p)))),
    ]

    fig = Figure()
    for (I, (v, what)) in pairs(collect(Iterators.product(views, whats)))
        pos = fig[Tuple(I)...]
        Axis(pos[1,1]; title="$(v.name), $(what.what|>typeof|>nameof) (θ=$(v.θ))", aspect=DataAspect(), width=300, height=300)
        jet = @set S.geometry_axis(jet0) = axis_for_viewing_angle(v.θ)
        img = render_field(jet; extent=3, ν=1, what.what)
        plt = heatmap!(img; what.kwargs(maximum(img))...)
        Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_thin$(suffix).png"), fig)

    # Thick jet with optional adaptive supersampling (CPU only)
    if is_cpu
        for adaptive_supersampling in (false, 4)
            fig = Figure()
            jet0 = @set jet0.emission.electrons.Ca_ordered = 9/jet0.emission.electrons.sinavg_a
            for (I, (v, what)) in pairs(collect(Iterators.product(views, whats)))
                pos = fig[Tuple(I)...]
                Axis(pos[1,1]; title="$(v.name), $(what.what|>typeof|>nameof) (θ=$(v.θ))", aspect=DataAspect(), width=300, height=300)
                jet = @set S.geometry_axis(jet0) = axis_for_viewing_angle(v.θ)
                img = render_field(jet; extent=3, ν=1, what.what, adaptive_supersampling)
                plt = heatmap!(img; what.kwargs(maximum(img))...)
                Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
                hidespines!()
                hidedecorations!()
            end
            resize_to_layout!()
            save(joinpath(outdir, adaptive_supersampling!==false ? "bk_jet_thick_ss$(suffix).png" : "bk_jet_thick$(suffix).png"), fig)
        end
    end
end

function bk_jet_thick_options_image(; render_field, log_path, suffix="")
    log_section!(log_path, "bk_jet_thick_options_image")

    φj = 2u"°"
    θ = 7u"°"
    axis = SVector(sin(θ), 0, cos(θ))

    Bscale = Profiles.Axial(S.PowerLaw(-1; val0=1, s0=1))
    Bconfigs = (
        (name="tangled", B=S.BFieldSpec(Bscale, Directions.Scalar(), b -> S.FullyTangled(b))),
        (name="poloidal", B=S.BFieldSpec(Bscale, Directions.Axial(), identity)),
        (name="toroidal", B=S.BFieldSpec(Bscale, Directions.Toroidal(), identity)),
        (name="helical ψ=5°", B=S.BFieldSpec(Bscale, Directions.HelicalAT(5u"°"), identity)),
        (name="helical ψ=85°", B=S.BFieldSpec(Bscale, Directions.HelicalAT(85u"°"), identity)),
        (name="mixture κ=1", B=S.BFieldSpec(Bscale, Directions.Axial(), b -> S.TangledOrderedMixture(b; kappa=1))),
    )

    models = (
        (name="isotropic", electrons=S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=0.1)),
        (name="aniso η=1/100", electrons=S.AnisotropicPowerLawElectrons(; p=2.3, η=0.01, Cj=1, Ca=0.1)),
        (name="aniso η=100", electrons=S.AnisotropicPowerLawElectrons(; p=2.3, η=100, Cj=1, Ca=0.1)),
    )

    fig = Figure()
    for (r, m) in enumerate(models), (c, bc) in enumerate(Bconfigs)
        pos = fig[r, c]
        Axis(pos[1, 1]; title="$(m.name), $(bc.name)", aspect=DataAspect(), width=300, height=300)

        unsupported = m.electrons isa S.AnisotropicPowerLawElectrons && bc.B.wrap !== identity
        if unsupported
            hidespines!()
            hidedecorations!()
            Label(pos[1, 1], "unsupported"; tellwidth=false, tellheight=false)
            continue
        end

        jet0 = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
            emission = S.SynchrotronEmission(;
                ne = Profiles.Axial(S.PowerLaw(-2; val0=1, s0=1)),
                bc.B,
                m.electrons,
            ),
        ) |> S.prepare_for_computations
        jet0 = @set jet0.emission.electrons.Ca_ordered = 9 / jet0.emission.electrons.sinavg_a

        img = render_field(jet0; extent=3, ν=1, what=S.Intensity())
        plt = heatmap!(img; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img)))
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_thick_options$(suffix).png"), fig)
    fig
end


function main(; render_field, log_path, suffix="", is_cpu=true)
    bk_jet_image(; render_field, log_path, suffix, is_cpu)
    bk_jet_thick_options_image(; render_field, log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
    logs = init_script(@__FILE__)
    main(; render_field=make_render_cpu(logs.log_cpu), log_path=logs.log_cpu)
    if HAS_METAL
        main(; render_field=make_render_metal(logs.log_metal), log_path=logs.log_metal, suffix="_metal", is_cpu=false)
    end
end

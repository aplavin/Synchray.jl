include("helpers.jl")


function conical_jet_polarization_evpa_image(; render_field, log_path, suffix="")
    log_section!(log_path, "conical_jet_polarization_evpa_image")

    φj = 2u"°"
    θ = 7u"°"
    axis = SVector(sin(θ), 0, cos(θ))

    Bscale = Profiles.Axial(S.PowerLaw(-1; val0=1, s0=1))
    Bconfigs = (
        (name="poloidal", B=S.BFieldSpec(Bscale, Directions.Axial(), identity)),
        (name="toroidal", B=S.BFieldSpec(Bscale, Directions.Toroidal(), identity)),
        (name="helical ψ=45°", B=S.BFieldSpec(Bscale, Directions.HelicalAT(45u"°"), identity)),
    )

    electrons = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=0.1)

    fig = Figure()

    # Row 1: Ordered fields
    for (c, bc) in enumerate(Bconfigs)
        pos = fig[1, c]
        Axis(pos[1, 1]; title="$(bc.name), ordered", aspect=DataAspect(), width=350, height=350)

        jet0 = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
            emission = S.SynchrotronEmission(;
                ne = Profiles.Axial(S.PowerLaw(-2; val0=1, s0=1)),
                bc.B,
                electrons,
            ),
        ) |> S.prepare_for_computations
        jet0 = @set jet0.emission.electrons.Ca_ordered = 9 / jet0.emission.electrons.sinavg_a

        img_IQU = render_field(jet0; extent=3, ν=1, what=S.IntensityIQU())
        img_I = getproperty.(img_IQU, :I)

        plt = heatmap!(img_I; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img_I)))
        evpa_ticks!(img_IQU; step=5, min_I_frac=1e-5)
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end

    # Row 2: TangledOrderedMixture fields (κ=2)
    κ = 2.0
    for (c, bc) in enumerate(Bconfigs)
        pos = fig[2, c]
        Axis(pos[1, 1]; title="$(bc.name), κ=$(κ)", aspect=DataAspect(), width=350, height=350)

        B_mixed = @set $(bc.B).wrap = b -> S.TangledOrderedMixture(b, κ)

        jet0 = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
            emission = S.SynchrotronEmission(;
                ne = Profiles.Axial(S.PowerLaw(-2; val0=1, s0=1)),
                B = B_mixed,
                electrons,
            ),
        ) |> S.prepare_for_computations
        jet0 = @set jet0.emission.electrons.Ca_ordered = 9 / jet0.emission.electrons.sinavg_a

        img_IQU = render_field(jet0; extent=3, ν=1, what=S.IntensityIQU())
        img_I = getproperty.(img_IQU, :I)

        plt = heatmap!(img_I; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img_I)))
        evpa_ticks!(img_IQU; step=5, min_I_frac=1e-5)
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end

    resize_to_layout!()
    save(joinpath(outdir, "conical_jet_polarization_evpa$(suffix).png"), fig)
    fig
end

function jet_combined_velocity_helicalrt_image(; render_field, log_path, suffix="")
    log_section!(log_path, "jet_combined_velocity_helicalrt_image")

    φj = 2u"°"
    θ = 7u"°"
    axis = SVector(sin(θ), 0, cos(θ))

    Bscale = Profiles.Axial(S.PowerLaw(-1; val0=1, s0=1))

    Bconfigs = (
        (name="HelicalRT ψ=30°", B=S.BFieldSpec(Bscale, Directions.HelicalRT(30u"°"), identity)),
        (name="HelicalRT ψ=60°", B=S.BFieldSpec(Bscale, Directions.HelicalRT(60u"°"), identity)),
    )

    velocity_configs = (
        (name="radial only",
         velocity=S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995))),
        (name="radial + rotation",
         velocity=S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.99)) +
                  S.VelocitySpec(Directions.Toroidal(), S.beta, Profiles.RigidRotation(; β_ref=0.02, ρ_ref=1.0))),
    )

    electrons = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=0.1)

    fig = Figure()
    for (r, vc) in enumerate(velocity_configs), (c, bc) in enumerate(Bconfigs)
        pos = fig[r, c]
        Axis(pos[1, 1]; title="$(vc.name), $(bc.name)", aspect=DataAspect(), width=350, height=350)

        jet = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            vc.velocity,
            emission = S.SynchrotronEmission(;
                ne = Profiles.AxialTransverse(S.PowerLaw(-2; val0=1, s0=1), η -> exp(-4 * η^2)),
                bc.B,
                electrons,
            ),
        ) |> S.prepare_for_computations

        img_IQU = render_field(jet; extent=3, ν=1, what=S.IntensityIQU())
        img_I = getproperty.(img_IQU, :I)

        plt = heatmap!(img_I; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img_I)))
        evpa_ticks!(img_IQU; step=5, min_I_frac=1e-5)
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "jet_combined_velocity_helicalrt$(suffix).png"), fig)
    fig
end


function main(; render_field, log_path, suffix="")
    conical_jet_polarization_evpa_image(; render_field, log_path, suffix)
    jet_combined_velocity_helicalrt_image(; render_field, log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
    logs = init_script(@__FILE__)
    main(; render_field=make_render_cpu(logs.log_cpu), log_path=logs.log_cpu)
    if HAS_METAL
        main(; render_field=make_render_metal(logs.log_metal), log_path=logs.log_metal, suffix="_metal")
    end
end

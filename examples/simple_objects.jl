include("helpers.jl")


function moving_ellipsoid_image(; render_field, log_path, suffix="")
    log_section!(log_path, "moving_ellipsoid_image")

    βs = [0.0 0.7; 0.9 0.99]
    fig = Figure()
    for (I, β) in pairs(βs)
        pos = fig[Tuple(I)...]
        Axis(pos[1,1]; title="β=$(β)", aspect=DataAspect(), width=300, height=300)
        plt = heatmap!(
            render_field(
                S.MovingUniformEllipsoid(
                    center=S.FourPosition(5.0, 0.0, 0.0, 5.0),
                    sizes=SVector(4, 1, 1),
                    u0=S.FourVelocity(SVector(β, 0, 0)),
                    jν=1,
                    αν=1e-5,
                );
                extent=2.0,
                ν=1.0,
                t=0.0,
                what=S.Intensity(),
            );
            colormap=:magma,
        )
        Colorbar(pos[1,2], plt)
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "moving_ellipsoid$(suffix).png"), fig)
    fig
end

function synchrotron_sphere_image(; render_field, log_path, suffix="")
    log_section!(log_path, "synchrotron_sphere_image")

    regimes = (
        (name="thin", Ca=1e-3),
        (name="intermediate", Ca=0.5),
        (name="thick", Ca=1e3),
    )

    ν = 1.0
    extent = 2.0

    make_obj(Ca) = S.UniformSynchrotronSphere(
        center=S.FourPosition(5.0, 0.0, 0.0, 5.0),
        radius=1,
        u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)),
        ne0=1.0,
        B0=1.0,
        electrons=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca),
    )

    fields = (
        (name="Intensity", what=S.Intensity(), colormap=:magma),
        (name="τ", what=S.OpticalDepth(), colormap=:viridis),
        (name="α", what=S.SpectralIndex(), colormap=:balance, colorrange=(-3, 3)),
    )

    fig = Figure()
    for (r, reg) in enumerate(regimes), (c, f) in enumerate(fields)
        pos = fig[r,c]
        Axis(pos[1,1]; title="$(reg.name): $(f.name)", aspect=DataAspect(), width=300, height=300)
        plt = heatmap!(render_field(make_obj(reg.Ca); extent, ν, what=f.what); f.colormap, colorrange=get(f, :colorrange, Makie.automatic))
        Colorbar(pos[1,2], plt)
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "synchrotron_sphere$(suffix).png"), fig)
    fig
end


function main(; render_field, log_path, suffix="")
    moving_ellipsoid_image(; render_field, log_path, suffix)
    synchrotron_sphere_image(; render_field, log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
    logs = init_script(@__FILE__)
    main(; render_field=make_render_cpu(logs.log_cpu), log_path=logs.log_cpu)
    if HAS_METAL
        main(; render_field=make_render_metal(logs.log_metal), log_path=logs.log_metal, suffix="_metal")
    end
end

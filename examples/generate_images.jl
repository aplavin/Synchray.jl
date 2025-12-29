import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

outdir = joinpath(@__DIR__, "images")
isdir(outdir) || mkpath(outdir)

const render_time_log_path = joinpath(@__DIR__, "render_times.txt")
open(render_time_log_path, "w") do _ end

using Synchray
import Synchray as S
using CairoMakie
using RectiGrids
using AxisKeysExtra
using PyFormattedStrings


function render_field(obj; extent, nz=256, npx=512, ν, t=0.0, what=S.Intensity())
    cam = S.CameraZ(; xys=grid(SVector, range(0±extent, npx), range(0±extent, npx)), nz, ν, t)

    result = open(render_time_log_path, "a") do io
        redirect_stdout(io) do
            S.render(cam, obj, what)
            label = f"{nameof(typeof(obj)):27s}, {nameof(typeof(what)):15s}, {npx}px × {nz}"
            @time label S.render(cam, obj, what)
        end
    end

    result
end

function moving_ellipsoid_image()
    βs = [0.0 0.7; 0.9 0.99]
    fig = Figure(; size=(420 * 2, 360 * 2))
    for (I, β) in pairs(βs)
        pos = fig[Tuple(I)...]
        Axis(pos[1,1]; title="β=$(β)", aspect=DataAspect())
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
    save(joinpath(outdir, "moving_ellipsoid.png"), fig)
    fig
end

function synchrotron_sphere_image()
    # Two optical-depth regimes by varying absorption normalization.
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
        model=S.PowerLawElectrons(; p=2.5, Cj=1.0, Ca),
    )

    fields = (
        (name="Intensity", what=S.Intensity(), colormap=:magma),
        (name="τ", what=S.OpticalDepth(), colormap=:viridis),
        (name="α", what=S.SpectralIndex(), colormap=:balance, colorrange=(-3, 3)),
    )

    fig = Figure(; size=(420 * length(fields), 360 * length(regimes)))
    for (r, reg) in enumerate(regimes), (c, f) in enumerate(fields)
        pos = fig[r,c]
        Axis(pos[1,1]; title="$(reg.name): $(f.name)", aspect=DataAspect())
        plt = heatmap!(render_field(make_obj(reg.Ca); extent, ν, what=f.what); f.colormap, colorrange=get(f, :colorrange, Makie.automatic))
        Colorbar(pos[1,2], plt)
        hidespines!()
        hidedecorations!()
    end
    save(joinpath(outdir, "synchrotron_sphere.png"), fig)
    fig
end

function main()
    moving_ellipsoid_image()
    synchrotron_sphere_image()
end

main()

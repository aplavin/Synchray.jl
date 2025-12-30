import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

outdir = joinpath(@__DIR__, "images")
isdir(outdir) || mkpath(outdir)

const render_time_log_path = joinpath(@__DIR__, "render_times.txt")
open(render_time_log_path, "w") do _ end

using Synchray
import Synchray as S
using MakieExtra; import CairoMakie
using RectiGrids
using AxisKeysExtra
using PyFormattedStrings
using Accessors


function render_field(obj; extent, nz=256, npx=256, ν, t=0.0, what=S.Intensity())
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
    save(joinpath(outdir, "synchrotron_sphere.png"), fig)
    fig
end

function bk_jet_image()
    φj = 4u"°"

    axis_for_viewing_angle(θ) = SVector(sin(θ), 0.0, cos(θ))

    jet0 = S.ConicalBKJet(;
        axis=SVector(NaN, NaN, NaN),
        φj,
        s=(1e-3 .. 50),
        s0=1.0,
        ne0=1.0,
        B0=1.0,
        speed_profile=(η -> (S.beta, 0.995)),
        model=S.PowerLawElectrons(; p=2.3, Cj=1.0, Ca=0.1),
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
        jet = @set jet0.axis = axis_for_viewing_angle(v.θ)
        img = render_field(jet; extent=3, ν=1, what.what)
        plt = heatmap!(img; what.kwargs(maximum(img))...)
        Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_thin.png"), fig)

    fig = Figure()
    jet0 = @set jet0.model.Ca = 9.0
    for (I, (v, what)) in pairs(collect(Iterators.product(views, whats)))
        pos = fig[Tuple(I)...]
        Axis(pos[1,1]; title="$(v.name), $(what.what|>typeof|>nameof) (θ=$(v.θ))", aspect=DataAspect(), width=300, height=300)
        jet = @set jet0.axis = axis_for_viewing_angle(v.θ)
        img = render_field(jet; extent=3, ν=1, what.what)
        plt = heatmap!(img; what.kwargs(maximum(img))...)
        Colorbar(pos[1,2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_thick.png"), fig)
end

function bk_jet_1_knot_snapshots_image()
    φj = 4u"°"
    θ = 3 * φj  # "small viewing angle", same as in bk_jet_image

    base_jet = S.ConicalBKJet(;
        axis = SVector(sin(θ), 0, cos(θ)),
        φj,
        s=(1e-3 .. 50),
        s0=1,
        ne0=1,
        B0=1,
        speed_profile=(η -> (S.beta, 0.995)),
        model=S.PowerLawElectrons(; p=2.3, Cj=1, Ca=1),
    )

    knot = let
        x_c0 = S.FourPosition(0, (0.1 * base_jet.axis)...)
        u = S.FourVelocity(0.995 * base_jet.axis)
        sizing = S.CrossSectionKnotSizing(0.3, 0.5)
        S.InertialEllipsoidalKnot(; x_c0, u, sizing, profile_ne=S.GaussianBump(100), profile_B=nothing)
    end

    jet = S.ConicalBKJetWithPatterns(base_jet, (knot,)) |> S.prepare_for_computations

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
    save(joinpath(outdir, "bk_jet_1_knot_snapshots.png"), fig)
    fig
end

function main()
    moving_ellipsoid_image()
    synchrotron_sphere_image()
    bk_jet_image()
    bk_jet_1_knot_snapshots_image()
end

main()

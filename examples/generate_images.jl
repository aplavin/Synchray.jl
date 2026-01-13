import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()

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


function evpa_ticks!(img_IQU; step=16, min_I_frac=0.03, color=:white)
    Imax = maximum(:I, img_IQU)

    subimg = @p let 
        CartesianIndices(img_IQU)
        first(__):CartesianIndex(step, step):last(__)
        img_IQU[__]
    end

    vals = @p let
        with_axiskeys(subimg)
        filter(((xy, s),) -> s.I ≥ min_I_frac * Imax)
    end

    scatter!(
        FPlot(vals, ((xy, s),) -> SVector(xy...); rotation=((xy, s),) -> S.evpa(s));
        marker=:hline, color)
end


function render_field(obj; extent, nz=256, npx=256, ν, t=0.0, what=S.Intensity(), adaptive_supersampling=false)
    cam = S.CameraZ(; xys=grid(SVector, x=range(0±extent, npx), y=range(0±extent, npx)), nz, ν, t)

    result = open(render_time_log_path, "a") do io
        redirect_stdout(io) do
            S.render(cam, obj, what; adaptive_supersampling)
            label = f"{nameof(typeof(obj)):27s}, {nameof(typeof(what)):15s}, {npx}px × {nz}"
            @time label S.render(cam, obj, what; adaptive_supersampling)
        end
    end

    result
end

function moving_ellipsoid_image()
    open(io -> write(io, "moving_ellipsoid_image()\n"), render_time_log_path, "a")

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
    open(io -> write(io, "synchrotron_sphere_image()\n"), render_time_log_path, "a")

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
        model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca),
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
    open(io -> write(io, "bk_jet_image()\n"), render_time_log_path, "a")

    φj = 4u"°"
    axis_for_viewing_angle(θ) = SVector(sin(θ), 0.0, cos(θ))

    jet0 = S.EmissionRegion(
        geometry = Geometries.Conical(; axis=SVector(NaN, NaN, NaN), φj, z=1e-3 .. 50),
        ne = Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
        B = S.BFieldSpec(Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), Directions.Scalar(), b->S.FullyTangled(b)),
        velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
        model = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1.0, Ca=0.1),
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
    save(joinpath(outdir, "bk_jet_thin.png"), fig)

    for adaptive_supersampling in (false, 4)
        fig = Figure()
        jet0 = @set jet0.model.Ca_ordered = 9/jet0.model.sinavg_a
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
        save(joinpath(outdir, adaptive_supersampling!==false ? "bk_jet_thick_ss.png" : "bk_jet_thick.png"), fig)
    end
end

function bk_jet_1_knot_snapshots_image()
    open(io -> write(io, "bk_jet_1_knot_snapshots_image()\n"), render_time_log_path, "a")

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
        ne = S.Profiles.Modified(S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)), knot),
        B = S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.beta, S.Profiles.Constant(0.995)),
        model = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=1),
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
    save(joinpath(outdir, "bk_jet_1_knot_snapshots.png"), fig)
    fig
end

function bk_jet_thick_options_image()
    open(io -> write(io, "bk_jet_thick_options_image()\n"), render_time_log_path, "a")

    φj = 2u"°"
    θ = 7u"°"
    axis = SVector(sin(θ), 0, cos(θ))

    Bscale = Profiles.Axial(S.PowerLaw(-1; val0=1, s0=1))
    Bconfigs = (
        (name="tangled", B=S.BFieldSpec(Bscale, Directions.Scalar(), b -> S.FullyTangled(b))),
        (name="poloidal", B=S.BFieldSpec(Bscale, Directions.Axial(), identity)),
        (name="toroidal", B=S.BFieldSpec(Bscale, Directions.Toroidal(), identity)),
        (name="helical ψ=5°", B=S.BFieldSpec(Bscale, Directions.Helical(5u"°"), identity)),
        (name="helical ψ=85°", B=S.BFieldSpec(Bscale, Directions.Helical(85u"°"), identity)),
        (name="mixture κ=1", B=S.BFieldSpec(Bscale, Directions.Axial(), b -> S.TangledOrderedMixture(b; kappa=1))),
    )

    models = (
        (name="isotropic", model=S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=0.1)),
        (name="aniso η=1/100", model=S.AnisotropicPowerLawElectrons(; p=2.3, η=0.01, Cj=1, Ca=0.1)),
        (name="aniso η=100", model=S.AnisotropicPowerLawElectrons(; p=2.3, η=100, Cj=1, Ca=0.1)),
    )

    fig = Figure()
    for (r, m) in enumerate(models), (c, bc) in enumerate(Bconfigs)
        pos = fig[r, c]
        Axis(pos[1, 1]; title="$(m.name), $(bc.name)", aspect=DataAspect(), width=300, height=300)

        unsupported = m.model isa S.AnisotropicPowerLawElectrons && bc.B.wrap !== identity
        if unsupported
            hidespines!()
            hidedecorations!()
            Label(pos[1, 1], "unsupported"; tellwidth=false, tellheight=false)
            continue
        end

        jet0 = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            ne = Profiles.Axial(S.PowerLaw(-2; val0=1, s0=1)),
            bc.B,
            velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
            m.model,
        ) |> S.prepare_for_computations
        jet0 = @set jet0.model.Ca_ordered = 9 / jet0.model.sinavg_a

        img = render_field(jet0; extent=3, ν=1, what=S.Intensity())
        plt = heatmap!(img; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img)))
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "bk_jet_thick_options.png"), fig)
    fig
end

function conical_jet_polarization_evpa_image()
    open(io -> write(io, "conical_jet_polarization_evpa_image()\n"), render_time_log_path, "a")

    φj = 2u"°"
    θ = 7u"°"
    axis = SVector(sin(θ), 0, cos(θ))

    Bscale = Profiles.Axial(S.PowerLaw(-1; val0=1, s0=1))
    Bconfigs = (
        (name="poloidal", B=S.BFieldSpec(Bscale, Directions.Axial(), identity)),
        (name="toroidal", B=S.BFieldSpec(Bscale, Directions.Toroidal(), identity)),
        (name="helical ψ=45°", B=S.BFieldSpec(Bscale, Directions.Helical(45u"°"), identity)),
    )

    model = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1, Ca=0.1)

    fig = Figure()
    for (c, bc) in enumerate(Bconfigs)
        pos = fig[1, c]
        Axis(pos[1, 1]; title="$(bc.name), I + EVPA", aspect=DataAspect(), width=350, height=350)

        jet0 = S.EmissionRegion(;
            geometry = Geometries.Conical(; axis, φj, z=1e-3 .. 50),
            ne = Profiles.Axial(S.PowerLaw(-2; val0=1, s0=1)),
            bc.B,
            velocity = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.995)),
            model,
        ) |> S.prepare_for_computations
        jet0 = @set jet0.model.Ca_ordered = 9 / jet0.model.sinavg_a

        img_IQU = render_field(jet0; extent=3, ν=1, what=S.IntensityIQU())
        img_I = getproperty.(img_IQU, :I)

        plt = heatmap!(img_I; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img_I)))
        evpa_ticks!(img_IQU; step=5, min_I_frac=1e-5)
        Colorbar(pos[1, 2], plt; tickformat=EngTicks(:symbol))
        hidespines!()
        hidedecorations!()
    end
    resize_to_layout!()
    save(joinpath(outdir, "conical_jet_polarization_evpa.png"), fig)
    fig
end

function main()
    moving_ellipsoid_image()
    synchrotron_sphere_image()
    bk_jet_image()
    bk_jet_1_knot_snapshots_image()
    bk_jet_thick_options_image()
    conical_jet_polarization_evpa_image()
end

main()

# using Pkg
# Pkg.activate(@__DIR__)
# Pkg.resolve()

using DataManipulation
using MakieExtra; import GLMakie
using Unitful
using AxisKeysExtra, RectiGrids
using StructArrays
using AccessibleModels; using AccessibleModels: P
using IntervarrrSets
using DistributionsExtra
using Interpolations
using Rotations
using PyFormattedStrings
using OhMyThreads: tmap
using Metal
using SQLCollections, SQLite

MakieExtra.show_gl_icon_in_dock()

using Synchray
import Synchray as S
import FieldLIC as LIC



function _SA_to_regular(sa::StructArray)
	map((xs...) -> eltype(sa)(Tuple(xs)), getproperties(sa)...)
end

_tmap(f, X::KeyedArray) = @set AxisKeys.keyless_unname(X) = tmap(f, X)


sstate = let
    db = SQLite.DB(joinpath(@__DIR__, "_interactive_states.db"))
    SQLDictionary{String, Float64}(db, :basic_nounits_state)
end


fig = Figure(size=(1600, 1000))

params, = SliderGrid(fig[1,1][1,1], AccessibleModel((;
    img=(
        npix=P(discreterange(log, 16..500, length=50)),
        nz=P(discreterange(log, 2..300, length=50)),
        dynrange=P(LogUniform(1..1e5), 1e3),
        ν1=P(LogUniform(0.1..10), 1),
    ),
    time=(
        t=P(-0.1..10),
        maxt=P(LogUniform(3..1e3)),
    ),
    geom=(
        z=0.1..1500,
        viewing_ang=P(LogUniform(0.1..80))u"°",
        opening_ang=P(LogUniform(0.1..80))u"°",
        ctrjet=P([false, true], false),
        β_φ=P(0..0.01, 0),
    ),
    B=(
        B0=P(LogUniform(1e-5..1e5)),
        helixψ=P(0..90)u"°",
        ordered=P(0..1, 0),
    ),
    electrons=(
        p=P(2:0.25:3, 2.5),
        ne0=P(LogUniform(1e-5..1e5)),
        anis=P([false, true], false),
        anis_η=P(LogUniform(1e-5..1e5)),
    ),
    knot=(
        enabled=P([false, true], false),
        wfrac=P(0..1, 0.01),
        elong=P(0..5, 1),
        γ=P(1.01..50),
        mul=P(LogUniform(5..1e3)),
    ),
    nozzle=(
        enabled=P([false, true], false),
        wfrac=P(0..1, 0.2),
        period=P(LogUniform(0.1..1e3)),
        mul=P(LogUniform(1..1e7))
    )
), AccessibleModels.Auto()); state=sstate, rowgap=2)
on(events(fig).tick) do tick
    params[].time.maxt < 3.1 && return
	params[] = @set $(params[]).time.t = mod(tick.time, 0..5)/5 * params[].time.maxt
end

γpoints = Observable([(x=0., γ=10.), (x=1., γ=10.)])
ip = InteractivePoints(γpoints)
axplot(scatter, widgets=[ip])(fig[1,1][2,1], (@lift FPlot($γpoints, (@o _.x), (@o _.γ))); axis=(;height=250, limits=(0..1, 0..50), title="Transverse velocity profile"))
γ_cross = @lift Profiles.LinearInterp(@p $γpoints map(_.x => _.γ) Tuple)::Any
lines!(0..1, @lift (x->$γ_cross(x))::Any)

camera = @lift S.CameraZ(;
    xys=grid(SVector,
            x=range(-5..30, length=2*$params.img.npix) |> LinRange,
            y=range(-5..5, length=$params.img.npix) |> LinRange),
    ν=$params.img.ν1, $params.img.nz, $params.time.t,
    mapfunc=_tmap)

jet = @lift let
    axis = RotY($params.geom.ctrjet ? π - $params.geom.viewing_ang : $params.geom.viewing_ang) * SVector(0,0,1)

    knot = S.Patterns.EllipsoidalKnot(
        x_c0=S.FourPosition(0.0, 0.0, 0.0, 0.0),
        u=construct(S.FourVelocity, S.gamma=>$params.knot.γ, S.direction=>axis),
        sizing=S.Patterns.CrossSectionSizing($params.knot.wfrac, $params.knot.elong),
        profile=S.Patterns.GaussianBump($params.knot.mul),
    )

    nozzle  = S.Patterns.PrecessingNozzle(
        θ_precession = $params.geom.opening_ang * (1 - $params.nozzle.wfrac),   # precession angle from axis
        θ_nozzle = $params.nozzle.wfrac * $params.geom.opening_ang,     # nozzle cone half-angle
        period = $params.nozzle.period,        # precession period in code time
        β_flow = S._beta_from_gamma($γ_cross(1)),      # flow speed
        profile = S.Patterns.complement(S.Patterns.TophatBump($params.nozzle.mul)),
    )


    ne_base = S.Profiles.Axial(S.PowerLaw(-2; val0=$params.electrons.ne0, s0=1.0))
    ne = $params.knot.enabled ? S.Profiles.Modified(ne_base, knot) :
         $params.nozzle.enabled ? S.Profiles.Modified(ne_base, nozzle) :
         ne_base

    velocity_along = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Transverse($γ_cross))

    region = S.EmissionRegion(;
        geometry = S.Geometries.Conical(; axis, φj=$params.geom.opening_ang, z=$params.geom.z),
        ne,
        B = S.BFieldSpec(
            S.Profiles.Axial(S.PowerLaw(-1; val0=$params.B.B0, s0=1.0)),
            $params.B.ordered > 0 ? S.Directions.HelicalAT($params.B.helixψ) : S.Directions.Scalar(),
            $params.B.ordered == 1 ? identity :
            $params.B.ordered == 0 ? b->S.FullyTangled(b) :
            b->S.TangledOrderedMixture(b; kappa=$params.B.ordered / (1 - $params.B.ordered)),
        ),
        velocity = $params.geom.β_φ == 0 ? velocity_along : velocity_along + S.VelocitySpec(S.Directions.Toroidal(), S.beta, S.Profiles.RigidRotation(β_ref=$params.geom.β_φ, ρ_ref=1.0)),
        model = $params.electrons.anis ?
            S.AnisotropicPowerLawElectrons(;$params.electrons.p, η=$params.electrons.anis_η, Cj=1.0, Ca=1.0) :
            S.IsotropicPowerLawElectrons(;$params.electrons.p, Cj=1.0, Ca=1.0),
    )

    S.prepare_for_computations(region)
end::Any

camera_gpu = @lift @p $camera S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, __.xys)
jet_gpu = @lift S.to_float_type(Float32, $jet)::Any

xy_sel = Observable(SVector(0., 0.))
yslice = @lift $xy_sel.y  # $params.img.yslice

Jcontrib = @lift let
    # x-range from camera FOV
    x_range = range(extrema(axiskeys($camera.xys, 1))..., 150)
    y = $yslice

    res = flatmap(x_range) do x_u
        x = x_u
        ray = S.RayZ(; x0=S.FourPosition($camera.t, x, y, 0.0), k=$camera.ν, nz=$camera.nz)
        prof = S.ray_contribution_profile_IQU($jet, ray)

        maxcontrib = maximum(s -> s.I, prof.dIν_to_obs; init=0.0)
        map(prof) do step
            x4 = step.x4
            r_lab = SVector(x4.x, x4.y, x4.z)
            r_local = S.rotate_lab_to_local($jet, r_lab)
            s_val = r_local[3]
            h_val = r_local[1]
            contrib = step.dIν_to_obs

            (;s=s_val, h=h_val, contrib=contrib.I / maxcontrib)
        end
    end |> StructArray
    sort(res; by=x->x.contrib)
end

let
    pos = fig[1:2,2][0,1]
    ax = Axis(pos[1,1];
        xlabel="s (along jet axis)",
        ylabel="h (perpendicular)",
        aspect=3,
        backgroundcolor=:black,)
    scatter!(ax,
        (@lift FPlot($Jcontrib, (@o _.s), (@o _.h), color=(@o _.contrib), markersize=10)),
        colormap=:inferno, colorrange=(0, 1))

    # Draw the (0,0) pixel's ray in jet-plane coordinates
    ray_path = @lift let
        xy = $xy_sel
        ray = S.RayZ(; x0=S.FourPosition($camera.t, xy..., 0.0), k=$camera.ν, nz=$camera.nz)
        axis_z = S.geometry_axis($jet).z
        z_range = 0±1e4
        S.ray_in_local_coords(ray, $jet; z_range)
    end
    lines!(ax, (@lift (@swiz $ray_path.zx)); color=:cyan, linewidth=1, linestyle=:dash, xautolimits=false, yautolimits=false)
end

img_iqu = @lift let
    (value, time) = @timed Array(S.render(camera_gpu[], $jet_gpu, S.IntensityIQU()))
    @info "" time
    KeyedArray(value; named_axiskeys(camera[].xys)...) |> StructArray
end
img_si = @lift KeyedArray(Array(S.render(camera_gpu[], $jet_gpu, S.SpectralIndex())); named_axiskeys(camera[].xys)...)

pos = fig[1:2,2][1,1]
ax, hm = image(pos[1,1],
    (@lift $img_iqu.I),
    colormap=:turbo,
    colorscale=(@lift SymLog(1/$params.img.dynrange * maximum($img_iqu.I))),
    colorrange=(@lift (0, 1) .* maximum($img_iqu.I)),
    interpolate=false,
    axis=(;title="Total intensity")
)
contour!((@lift $img_iqu.I), color=(:gray, 0.5), linewidth=1, levels=@lift @p maximum($img_iqu.I) maprange(log, (0.05/$params.img.dynrange * __)..__, length=30))
hlines!((@lift $yslice); color=(:white, 0.5), linestyle=:dash)
scatter!((@lift $xy_sel); color=:transparent, markersize=15, marker=:circle, strokewidth=1, strokecolor=:cyan)
Colorbar(pos[1,2], hm)

on(mouse_position_obs(ax; key=Mouse.left)) do pos
    xy_sel[] = SVector(pos...)
end

lines(fig[2,1][1,1], (@lift @p $img_iqu.I(x=Near(20.0))  __ ./ maximum(__)))

# pos = fig[1:2,2][2,1]
# ax, hm = image(pos[1,1],
#     (@lift LIC.polarization_lic($img_iqu |> _SA_to_regular |> permutedims;
#         # combiner=LIC.BlendLinesOverlay(; lic_color=Makie.RGB(1.0, 1.0, 1.0),
#         #        colormap=Makie.ColorSchemes.turbo,
#         # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
#         # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
#         combiner=LIC.DraperyOverlay(; alpha=0.8, dark=false,
#             colormap=Makie.ColorSchemes.turbo,
#             colorscale=identity,
#             colorrange=(0, 0.7)),
#         # combiner=LIC.DraperyOverlay(; alpha=0.8, dark=true,
#         #        colormap=Makie.ColorSchemes.turbo,
#         # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
#         # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
#         # combiner=LIC.MultiplicativeDrapery(; alpha_dark=0.5, alpha_light=1,
#         #        colormap=Makie.ColorSchemes.turbo,
#         # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
#         # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
#         scale=2,
#         noise_generator=LIC.SparseNoise(density=0.02, seed=42),
#         # noise_generator=LIC.UpscaledNoise(scale=3, seed=42),
#         lic=let L=250
#             LIC.FLIC(float(L), ones(2L+1)/(2L+1); nematic=true)
#         end,
#         flux_threshold=1e-5,
#         max_pol_fraction=0.4,
#         pol_fraction_gamma=0.5,
#         base_quantity=:fraction
#     ) |> permutedims),
#     interpolate=false,
#     axis=(;title="Polarized fraction + orientation")
# )

pos = fig[1:2,2][3,1]
ax,plt = image(pos[1,1], img_si,
    colormap=:turbo, colorrange=(-1, 2.6),
    interpolate=false,
    axis=(;title="Spectral index")
)
Colorbar(pos[1,2], plt)

fig |> display

# @profview for _ in 1:10
#     @time params[] = params[]
# end

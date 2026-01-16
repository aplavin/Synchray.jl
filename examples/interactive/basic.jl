# using Pkg
# Pkg.activate(@__DIR__)
# Pkg.resolve()

using DataManipulation
using MakieExtra; import GLMakie
using Unitful, UnitfulAstro
using AxisKeysExtra, RectiGrids
using StructArrays
using AccessibleModels; using AccessibleModels: P
using IntervarrrSets
using DistributionsExtra
using Interpolations
using Rotations
using PyFormattedStrings
using OhMyThreads: tmap
using SQLCollections, SQLite

MakieExtra.show_gl_icon_in_dock()

using Synchray
import Synchray as S
import FieldLIC as LIC





@eval LIC begin
function combine(oc::DraperyOverlay, intensity, lic)
    int_interp = _interpolate_to_lic(intensity, lic)
    lo, hi = oc.colorrange
	llo, lhi = extrema(lic)
	_combine(oc::DraperyOverlay, int_interp, lic, lo, hi, llo, lhi; oc.colorscale, oc.colormap, oc.alpha, oc.dark)
end
function _combine(oc::DraperyOverlay, int_interp, lic, lo, hi, llo, lhi; colorscale, colormap, alpha, dark)
    map(int_interp, lic) do int, l
		l = (l - llo) / (lhi - llo)
        # Get base color from intensity
        int_clamped = clamp(int, lo, hi)
        t = (colorscale(int_clamped) - colorscale(lo)) / (colorscale(hi) - colorscale(lo))
        t = clamp(t, 0, 1)
        base_color = get(colormap, t)

        # LIC texture as grayscale overlay
        l_clamped = clamp(l, 0, 1)
        if dark
            # Dark drapery: blend toward black
            overlay_color = RGB(0.0, 0.0, 0.0)
        else
            # Light drapery: blend toward white
            overlay_color = RGB(1.0, 1.0, 1.0)
        end
		drapery_strength = l_clamped * alpha
		t == 0 && (drapery_strength *= 0)

        (1 - drapery_strength) * base_color + drapery_strength * overlay_color
    end
end
function combine(oc::BlendLinesOverlay, intensity, lic)
    int_interp = _interpolate_to_lic(intensity, lic)
    lo, hi = oc.colorrange
    map(int_interp, lic) do int, l
        # Map intensity to color via colormap
        int_clamped = clamp(int, lo, hi)
        t = (oc.colorscale(int_clamped) - oc.colorscale(lo)) / (oc.colorscale(hi) - oc.colorscale(lo))
        t = clamp(t, 0, 1)
        base_color = get(oc.colormap, t)

        # Blend with line color based on LIC
        l_clamped = t == 0 ? zero(l) : clamp(l, 0, 1)
        (1 - l_clamped) * base_color + l_clamped * oc.lic_color
    end
end
end




function _SA_to_regular(sa::StructArray)
	map((xs...) -> eltype(sa)(Tuple(xs)), getproperties(sa)...)
end

_tmap(f, X::KeyedArray) = @set AxisKeys.keyless_unname(X) = tmap(f, X)


sstate = let
    db = SQLite.DB(joinpath(@__DIR__, "_interactive_states.db"))
    SQLDictionary{String, Float64}(db, :basic_state)
end


fig = Figure(size=(1600, 1000))

params, = SliderGrid(fig[1,1][1,1], AccessibleModel((;
    img=(
        npix=P(discreterange(log, 16..500, length=50)),
        nz=P(discreterange(log, 2..300, length=50)),
        dynrange=P(LogUniform(1..1e5), 1e3),
        evpastep=P(1:20, 3),
        t=P(-0.1..10)u"yr",
        maxt=P(LogUniform(3..1e3))u"yr",
        ν1=P(LogUniform(1..1e3), 15)*u"GHz",
    ),
    geom=(
        z=(0.1..500)u"pc",
        # z=(P(LogUniform(1e-3..1e4), 1e-3)..P(LogUniform(1e-3..1e4), 1e3))u"pc",
        viewing_ang=P(LogUniform(0.1..80))u"°",
        opening_ang=P(LogUniform(0.1..80))u"°",
        ctrjet=P([false, true], false),
    ),
    B=(
        B0=P(LogUniform(1e-5..1e5))*u"Gauss",
        helixψ=P(0..90)u"°",
        ordered=P([false,true], false),
    ),
    electrons=(
        p=P(2:0.25:3, 2.5),
        ne0=P(LogUniform(1e-5..1e5))*u"cm^-3",
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
        period=P(LogUniform(0.1..1e3))u"yr",
        mul=P(LogUniform(1..1e7))
    )
), AccessibleModels.Auto()); state=sstate, rowgap=2)
on(events(fig).tick) do tick
    params[].img.maxt < 3.1u"yr" && return
	params[] = @set $(params[]).img.t = mod(tick.time, 0..5)/5 * params[].img.maxt
end

γpoints = Observable([(x=0., γ=10.), (x=1., γ=10.)])
ip = InteractivePoints(γpoints)
axplot(scatter, widgets=[ip])(fig[1,1][2,1], (@lift FPlot($γpoints, (@o _.x), (@o _.γ))); axis=(;height=250, limits=(0..1, 0..50), title="Transverse velocity profile"))
γ_cross = @lift Profiles.LinearInterp(@p $γpoints map(_.x => _.γ) Tuple)
lines!(0..1, @lift x->$γ_cross(x))

camera = @lift S.CameraZ(;
    xys=grid(SVector,
            x=range(-5..30, length=2*$params.img.npix) .* u"pc" |> LinRange,
            y=range(-5..5, length=$params.img.npix) .* u"pc" |> LinRange),
    ν=$params.img.ν1, $params.img.nz, $params.img.t,
    mapfunc=_tmap) |> ustrip

jet = @lift let
    axis = RotY($params.geom.ctrjet ? π - $params.geom.viewing_ang : $params.geom.viewing_ang) * SVector(0,0,1)
    
    knot = S.Patterns.EllipsoidalKnot(
        x_c0=S.FourPosition(0.0u"pc", 0.0u"pc", 0u"pc", 0u"pc"),
        # x_c0=S.FourPosition(0.0u"yr", 0.0u"pc", 0u"pc", 2u"pc"),
        u=construct(S.FourVelocity, S.gamma=>$params.knot.γ, S.direction=>axis),
        sizing=S.Patterns.CrossSectionSizing($params.knot.wfrac, $params.knot.elong),
        profile=S.Patterns.GaussianBump($params.knot.mul),
    )

    nozzle  = S.Patterns.PrecessingNozzle(
        θ_precession = $params.geom.opening_ang * (1 - $params.nozzle.wfrac),   # precession angle from axis
        θ_nozzle = $params.nozzle.wfrac * $params.geom.opening_ang,     # nozzle cone half-angle
        period = $params.nozzle.period,        # precession period in code time
        β_flow = √(1 - $γ_cross(1)^-2),      # flow speed
        profile = S.Patterns.TophatBump($params.nozzle.mul),
    )


    ne_base = S.Profiles.Axial(S.PowerLaw(-2; val0=$params.electrons.ne0, s0=1u"pc"))
    ne = $params.knot.enabled ? S.Profiles.Modified(ne_base, knot) :
         $params.nozzle.enabled ? S.Profiles.Modified(ne_base, nozzle) :
         ne_base
    region = S.EmissionRegion(;
        geometry = S.Geometries.Conical(; axis, φj=$params.geom.opening_ang, z=$params.geom.z),
        ne,
        B = S.BFieldSpec(
            S.Profiles.Axial(S.PowerLaw(-1; val0=$params.B.B0, s0=1u"pc")),
            $params.B.ordered ? S.Directions.Helical($params.B.helixψ) : S.Directions.Scalar(),
            $params.B.ordered ? identity : b->S.FullyTangled(b),
        ),
        velocity = S.VelocitySpec(S.Directions.Axial(), S.gamma, S.Profiles.Transverse($γ_cross)),
        model = $params.electrons.anis ?
            S.AnisotropicPowerLawElectrons(;$params.electrons.p, η=$params.electrons.anis_η) :
            S.IsotropicPowerLawElectrons(;$params.electrons.p),
    )
    
    S.prepare_for_computations(region |> ustrip)
end::Any

JG = @lift @p grid(SVector, x=range(0±(500*tan($jet.geometry.φj)), 3*$params.img.npix)u"pc", y=[0]u"pc", z=range(0..500, 3*$params.img.npix)u"pc") dropdims(dims=:y) permutedims() @modify(x->ustrip.(x), axiskeys(__)[∗])
Jimg = @lift map($JG) do rj
    r = S.rotate_local_to_lab($jet, S._u_to_code(rj, S.UCTX.L0))
    x = S.event_on_camera_ray($camera, r)
    (ne=S.is_inside($jet, x) * S.electron_density($jet, x),)
end |> StructArray
# Jimg[] |> display
# camera[]

let img = @lift $Jimg.ne
    pos = fig[1:2,2][0,1]
    image(pos[1,1], img, colormap=:turbo, colorscale=(@lift SymLog(1/$params.img.dynrange * maximum($img))), 
        axis=(;title="Electron density (jet crossection)", aspect=3))
end
# lines!((@lift sightline($jet, $params.z)); to_xy_attrs(autolimits=false)...)

# # let img = @lift map(xj -> S.jet_at(Val(:j_nu_contrib), $jet, xj, params[].ν_obs, params[].z) |> ustrip, JG[])
# # 	image(fig[2,0][1,1], img, colormap=:turbo, colorscale=(@lift SymLog(1e-8maximum($img))))
# # end

img_iqu = @lift let
    (value, time) = @timed S.withunits(S.render, camera[], $jet, S.IntensityIQU())
    @info "" time
    map(v -> ustrip.(u"Jy/mas^2", v), value) |> StructArray
end
img_si = @lift S.render($camera, $jet, S.SpectralIndex())

# 	img = @lift first.($img_i_si)
# 	img_si = @lift last.($img_i_si) .|> NoUnits

    # lines(fig[2,0], @lift $img_all.:1(y=Near(0u"pc")))

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
    # evpa_ticks!(img_iqu; step=(@lift $params.img.evpastep), min_I_frac=1e-5)
    Colorbar(pos[1,2], hm)

lines(fig[2,1][1,1], (@lift @p $img_iqu.I(x=Near(20u"pc"))  __ ./ maximum(__)))

    pos = fig[1:2,2][2,1]
    ax, hm = image(pos[1,1], 
        (@lift LIC.polarization_lic($img_iqu |> _SA_to_regular |> permutedims;
            # combiner=LIC.BlendLinesOverlay(; lic_color=Makie.RGB(1.0, 1.0, 1.0),
            #        colormap=Makie.ColorSchemes.turbo,
            # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
            # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
            combiner=LIC.DraperyOverlay(; alpha=0.8, dark=false,
                colormap=Makie.ColorSchemes.turbo,
                colorscale=identity,
                colorrange=(0, 0.7)),
            # combiner=LIC.DraperyOverlay(; alpha=0.8, dark=true,
            #        colormap=Makie.ColorSchemes.turbo,
            # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
            # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
            # combiner=LIC.MultiplicativeDrapery(; alpha_dark=0.5, alpha_light=1,
            #        colormap=Makie.ColorSchemes.turbo,
            # 	colorscale=SymLog(1/$params.img.dynrange * maximum($img_iqu.I)),
            # 	colorrange=(0, 1) .* maximum($img_iqu.I),),
            scale=2,
            noise_generator=LIC.SparseNoise(density=0.02, seed=42),
            # noise_generator=LIC.UpscaledNoise(scale=3, seed=42),
            lic=let L=250
                LIC.FLIC(float(L), ones(2L+1)/(2L+1); nematic=true)
            end,
            flux_threshold=1e-5,
            max_pol_fraction=0.4,
            pol_fraction_gamma=0.5,
            base_quantity=:fraction
        ) |> permutedims),
        interpolate=false,
        axis=(;title="Polarized fraction + orientation")
    )

# 	on(mouse_position_obs(ax)) do pos
# 		xy_obs = SVector(pos[1], pos[2])*u"pc"
# 		(;xyz_jet) = closest_to_cone_axis_on_los(jet[], xy_obs)
# 		@reset xyz_jet.z /= 150
# 		# @info "" xy_obs xyz_jet
# 		knots[] = @set $(knots[])[1].x_inj = xyz_jet
# 	end

    pos = fig[1:2,2][3,1]
    ax,plt = image(pos[1,1], img_si,
        colormap=:turbo, colorrange=(-1, 2.6),
        interpolate=false,
        axis=(;title="Spectral index")
    )
    Colorbar(pos[1,2], plt)

fig |> display

# @time params[] = @set $(params[]).img.ν1 = params[].img.ν1;
# # @profview @time params[] = @set $(params[]).img.ν1 = params[].img.ν1;
# @profview for _ in 1:20
# 	@time params[] = @set $(params[]).img.ν1 = params[].img.ν1
# end

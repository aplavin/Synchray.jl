import Pkg
Pkg.activate(@__DIR__)

using Synchray
import Synchray as S
using MakieExtra; import CairoMakie
using RectiGrids
using AxisKeysExtra

# ===== Jet parameters =====
φj = 3u"°"
viewing_angles = [2u"°", 6u"°", 10u"°"]  # columns
period = 6.0  # precession period
γ = 10.0  # Lorentz factor

# ===== Create jets =====
println("Creating jets...")
jets = map(viewing_angles) do θ
    axis = SVector(sin(θ), 0.0, cos(θ))
    S.EmissionRegion(
        geometry = S.Geometries.Conical(; axis, φj, z = 0.1 .. 200.0),
        ne = S.Profiles.Modified(
            S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
            S.Patterns.PrecessingNozzle(
                θ_precession = φj * (1 - 1/5),   # precession angle from axis
                θ_nozzle = φj/5,     # nozzle cone half-angle
                period = period,     # precession period in code time
                β_flow = S._beta_from_gamma(γ),  # flow speed from Lorentz factor
                profile = S.Patterns.complement(S.Patterns.TophatBump(1e5))
            )
        ),
        B = S.BFieldSpec(
            S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)),
            S.Directions.Toroidal(),
            identity
        ),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Constant(γ)),
        model = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=10),
    ) |> S.prepare_for_computations
end

# ===== Video parameters =====
extent = 3.0
npx = 256
EVPA_STEP = 2
nz = 128
n_frames = 60
ν = 1.0
t_range = period * range(0, 1, length=n_frames+1)[1:end-1]

using OhMyThreads: tmap
_tmap(f, X::KeyedArray) = @set AxisKeys.keyless_unname(X) = tmap(f, X)

# ===== Helper function for EVPA ticks =====
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

# ===== Create figure with 3 columns =====
t = Observable(0.0)
fig = Figure(; size=(1800, 1400), backgroundcolor = :black)

println("Rendering initial frames...")

# Create panels for each viewing angle
for (i_col, (θ, jet)) in enumerate(zip(viewing_angles, jets))
    # Camera and rendering observables
    cam = @lift S.CameraZ(; xys=grid(SVector, x=range(2±extent, npx), y=range(0±extent, npx)), nz, ν, t=$t, mapfunc=_tmap)
    img_IQU = @lift S.render($cam, jet, S.IntensityIQU())
    img_I = @lift map(s -> s.I, $img_IQU)

    # Compute polarization observables
    pol_frac = @lift map(s -> sqrt(s.Q^2 + s.U^2) / (s.I + 1e-30), $img_IQU)

    # ===== Row 1: Total Intensity =====
    ax1 = Axis(fig[1, i_col];
        aspect = DataAspect(),
        title = "Total Intensity (θ=$(Int(round(rad2deg(θ))))°)",
        titlecolor = :gray,
        bottomspinecolor=:gray, topspinecolor=:gray, leftspinecolor=:gray, rightspinecolor=:gray
    )
    hidedecorations!(ax1)

    # Render initial frame to get color limits
    colorrange_I = maximum(img_I[]) .* (0, 1)
    colorscale_I = SymLog(3e-4*maximum(img_I[]))

    image!(ax1, img_I; colorscale=colorscale_I, colorrange=colorrange_I, colormap=:afmhot)

    # ===== Row 2: Polarization Fraction + EVPA =====
    ax2 = Axis(fig[2, i_col];
        aspect = DataAspect(),
        title = "Polarization Fraction + EVPA",
        titlecolor = :gray,
        bottomspinecolor=:gray, topspinecolor=:gray, leftspinecolor=:gray, rightspinecolor=:gray
    )
    hidedecorations!(ax2)

    plt = image!(ax2, pol_frac; colormap=:viridis, colorrange=(0, 0.7))
    if i_col == length(viewing_angles)
        Colorbar(fig[2, i_col+1], plt)
    end

    # EVPA ticks that update with time
    evpa_tick_obs = @lift begin
        img = $img_IQU
        Imax = maximum(:I, img)

        subimg = @p let
            CartesianIndices(img)
            first(__):CartesianIndex(EVPA_STEP, EVPA_STEP):last(__)
            img[__]
        end

        vals = @p let
            with_axiskeys(subimg)
            filter(((xy, s),) -> s.I ≥ 1e-4 * Imax)
        end

        vals
    end

    scatter!(ax2,
        (@lift FPlot($evpa_tick_obs, ((xy, s),) -> SVector(xy...); rotation=((xy, s),) -> S.evpa(s)));
        marker=:hline, color=:black)

    # ===== Row 3: Timeseries =====
    ax3 = Axis(fig[3, i_col];
        xlabel = "Time",
        ylabel = i_col == 1 ? "Flux" : "",
        xautolimitmargin=(0,0),
        bottomspinecolor=:gray, topspinecolor=:gray, leftspinecolor=:gray, rightspinecolor=:gray,
        backgroundcolor = :transparent
    )

    # Create secondary y-axis for EVPA
    ax3_evpa = Axis(fig[3, i_col];
        ylabel = i_col == length(viewing_angles) ? "EVPA" : "",
        yaxisposition = :right,
        yticks=Makie.AngularTicks(rad2deg(1), "°"),
        xautolimitmargin=(0,0),
        backgroundcolor = :transparent,
        ylabelcolor=:cyan, yticklabelcolor=:cyan,
        bottomspinecolor=:transparent, topspinecolor=:transparent,
        leftspinecolor=:transparent, rightspinecolor=:cyan
    )
    hidexdecorations!(ax3_evpa)

    # Current point observables (computed from current img_IQU)
    total_flux = @lift sum(s -> s.I, $img_IQU)
    total_pol_flux = @lift sum(s -> sqrt(s.Q^2 + s.U^2), $img_IQU)
    total_evpa = @lift begin
        sum_Q = sum(s -> s.Q, $img_IQU)
        sum_U = sum(s -> s.U, $img_IQU)
        mod(0.5 * atan(sum_U, sum_Q) + π/2, π) - π/2
    end
    timedata = obsmap(t, t_range, @lift (
        t=$t,
        flux=$total_flux,
        pol_flux=$total_pol_flux,
        pol_frac=$total_pol_flux / $total_flux,
        evpa=$total_evpa,
    ))
    fluxmul = 1 / maximum(x->x.flux, timedata) * 0.72
    tfplt = FPlot(timedata, (@o _.t))

    # Plot full timeseries using obsmap
    lines!(ax3, (@insert tfplt[2] = @o _.flux * fluxmul); color=:orange, linewidth=2, label="Total Flux (norm)")
    lines!(ax3, (@insert tfplt[2] = @o _.pol_frac); color=:magenta, linewidth=2, label="Pol Frac")
    lines!(ax3_evpa, (@insert tfplt[2] = @o _.evpa); color=:cyan, linewidth=2)

    # Current point markers
    scatter!(ax3, (@lift ($t, $total_flux * fluxmul)); color=:orange, markersize=15, strokewidth=2, strokecolor=:white)
    scatter!(ax3, (@lift ($t, $total_pol_flux / $total_flux)); color=:magenta, markersize=15, strokewidth=2, strokecolor=:white)
    scatter!(ax3_evpa, (@lift ($t, $total_evpa)); color=:cyan, markersize=15, strokewidth=2, strokecolor=:white)

    if i_col == 1
        axislegend(ax3; position=:lt, labelcolor=:white, backgroundcolor=(:black, 0.7))
    end
end

# Record video
println("Recording video...")
rec = Record(t, t_range; framerate=30)
save("single_jet_polarization.mp4", rec)

println("Video saved to single_jet_polarization.mp4")

# Create looped version
println("Creating looped version...")
run(`ffmpeg -y -stream_loop 4 -i single_jet_polarization.mp4 -c copy single_jet_polarization_loop.mp4`)
println("Looped video saved to single_jet_polarization_loop.mp4")
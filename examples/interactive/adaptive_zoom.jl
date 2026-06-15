# Adaptive zoom demo: pixel grid always covers the visible axis range at fixed resolution.
# Scroll to zoom in/out — the image re-renders at 256px in the new visible range.
# Two panels: ν=1 (left) and ν=3 (right), paraboloidal jet with adjustable shape.

using MakieExtra; import GLMakie
using RectiGrids
using AxisKeysExtra
using Rotations
using Metal
using DataManipulation
using Unitful
using AccessibleModels; using AccessibleModels: P
using DistributionsExtra

using Synchray
import Synchray as S

MakieExtra.show_gl_icon_in_dock()

# --- Layout ---

const NPIX = 256
fig = Figure(size=(1400, 500))

# --- Sliders ---

params, = SliderGrid(fig[0, 1:2], AccessibleModel((;
    a=P(0.05..1.0, 0.5),
    R_ref=P(LogUniform(1e-3..1e2), 1),
), AccessibleModels.Auto()); tellwidth=false, rowgap=2)

# --- Jet model (reactive) ---

jet = @lift let
    axis = RotY(1.5u"°") * SVector(0, 0, 1)

    j = S.EmissionRegion(;
        geometry = S.Geometries.Parabolic(; axis, R_ref=$params.R_ref, z_ref=50.0, a=$params.a, z=0.1..1000.0),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Transverse(Returns(5.0))),
        emission = S.SynchrotronEmission(;
            ne = S.Profiles.Axial(S.PowerLaw(-2*$params.a; val0=1.0, s0=1.0)),
            B = S.BFieldSpec(
                S.Profiles.Axial(S.PowerLaw(-$params.a; val0=1.0, s0=1.0)),
                S.Directions.Scalar(),
                b -> S.FullyTangled(b),
            ),
            electrons = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=10.0),
        ),
    )
    # j |> S.prepare_for_computations(
end::Any

jet_gpu = @lift S.to_float_type(Float32, $jet)::Any

# --- View state ---

view_rect = Observable(Rect2f(-5, -5, 35, 10))  # x: -5..30, y: -5..5

# --- Cameras at ν=1 and ν=3 ---

function make_camera(view_rect, ν)
    @lift let
        xi, yi = intervals($view_rect)
        S.CameraZ(;
            xys = grid(SVector,
                x = range(xi, length=2*NPIX) |> LinRange,
                y = range(yi, length=NPIX) |> LinRange),
            ν = Float32(ν), nz = 100, t = 0.0f0)
    end
end

camera1 = make_camera(view_rect, 1.0)
camera3 = make_camera(view_rect, 3.0)

camera1_gpu = @lift @p $camera1 S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))
camera3_gpu = @lift @p $camera3 S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))

# --- Render ---

img1 = @lift begin
    (value, elapsed) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera1_gpu, $jet_gpu, S.Intensity())))
    @info "rendered ν=1" elapsed
    value
end

img3 = @lift begin
    (value, elapsed) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera3_gpu, $jet_gpu, S.Intensity())))
    @info "rendered ν=3" elapsed
    value
end

# --- Plot panels ---

for (col, img, ν_label) in [(1, img1, "1"), (2, img3, "3")]
    ax, hm = image(fig[1, col], img;
        colormap = :afmhot,
        colorscale = (@lift let peak = maximum($img); SymLog(1e-1 * peak) end),
        colorrange = (@lift (0.0, maximum($img))),
        interpolate = false,
        axis = (;
            backgroundcolor=:black,
            xlabel = "x (code units)",
            ylabel = "y (code units)",
            title  = "ν = $ν_label",
        ),
    )
    Colorbar(fig[1, col][1, 2], hm)
end

# Sync zoom: bind ax.finallimits → view_rect (throttled)
connect!(view_rect, Makie.Observables.throttle(0.05, contents(fig[1, 1])[1].finallimits))

linkaxes!(contents(fig[1, 1])[1], contents(fig[1, 2])[1])

display(fig)

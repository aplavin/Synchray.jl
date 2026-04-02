# Adaptive zoom demo: pixel grid always covers the visible axis range at fixed resolution.
# Scroll to zoom in/out — the image re-renders at 256px in the new visible range.

using MakieExtra; import GLMakie
using RectiGrids
using AxisKeysExtra
using Rotations
using Metal
using DataManipulation
using Unitful

using Synchray
import Synchray as S

MakieExtra.show_gl_icon_in_dock()

# --- Fixed jet model ---

jet = let
    axis = RotY(1.5u"°") * SVector(0, 0, 1)

    S.prepare_for_computations(S.EmissionRegion(;
        geometry = S.Geometries.Conical(; axis, φj=1u"°", z=0.1..1000.0),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Transverse(Returns(5.0))),
        emission = S.SynchrotronEmission(;
            ne = S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
            B = S.BFieldSpec(
                S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)),
                S.Directions.Scalar(),
                b -> S.FullyTangled(b),
            ),
            electrons = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1.0),
        ),
    ))
end
jet_gpu = S.to_float_type(Float32, jet)

# --- View state: one Rect2f matching ax.finallimits ---

const NPIX = 256
view_rect = Observable(Rect2f(-5, -5, 35, 10))  # x: -5..30, y: -5..5

# --- Camera and render, both derived from view_rect ---

camera = @lift let
    xi, yi = intervals($view_rect)
    S.CameraZ(;
        xys = grid(SVector,
            x = range(xi, length=2*NPIX) |> LinRange,
            y = range(yi, length=NPIX) |> LinRange),
        ν = 1.0f0, nz = 100, t = 0.0f0)
end

camera_gpu = @lift @p $camera S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))

img = @lift begin
    (value, elapsed) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera_gpu, jet_gpu, S.Intensity())))
    @info "rendered" elapsed rect=view_rect[]
    value
end

# --- Plot (image() uses axis keys from the KeyedArray for physical coordinates) ---

fig = Figure(size=(900, 400))
ax, hm = image(fig[1, 1], img;
    colormap = :afmhot,
    colorrange = (@lift (0.0, maximum($img))),
    axis = (;
        backgroundcolor=:black,
        xlabel = "x (code units)",
        ylabel = "y (code units)",
        title  = "Adaptive zoom — scroll to zoom ($(NPIX) px always covers visible range)",
    ),
)
Colorbar(fig[1, 2], hm)

# --- Bind ax.finallimits → view_rect with throttle ---
# During zoom, Makie stretches the current image as preview.
# After 300 ms of no new events, view_rect updates → camera + render fire.

# connect!(view_rect, ax.finallimits)
connect!(view_rect, Makie.Observables.throttle(0.05, ax.finallimits))

display(fig)

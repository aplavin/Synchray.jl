import Pkg
Pkg.activate(@__DIR__)

using Synchray
import Synchray as S
using MakieExtra; import CairoMakie
using RectiGrids
using AxisKeysExtra

# ===== Setup EmissionRegion with precessing nozzle pattern =====
φj = 4u"°"
θ = 3 * φj  # viewing angle
axis = SVector(sin(θ), 0.0, cos(θ))

jet = S.EmissionRegion(
    geometry = S.Geometries.Conical(; axis, φj, z = 1.0 .. 50.0),
    ne = S.Profiles.Modified(
        S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
        S.Patterns.PrecessingNozzle(
            θ_precession = φj,   # precession angle from axis
            θ_nozzle = φj/5,     # nozzle cone half-angle
            period = 5.0,        # precession period in code time
            β_flow = 0.995,      # flow speed
            profile = S.Patterns.TophatBump(1e4)
        )
    ),
    B = S.BFieldSpec(
        S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)), 
        S.Directions.Helical(75u"°"), 
        identity
    ),
    velocity = S.VelocitySpec(S.Directions.Radial(), S.beta, S.Profiles.Constant(0.995)),
    model = S.IsotropicPowerLawElectrons(; p=2.3, Cj=1.0, Ca=0.1),
) |> S.prepare_for_computations


# ===== Render at three different times =====
extent = 3.0
npx = 128
nz = 256
ν = 1.0

# Three phases: 0, 1/3, 2/3 of the precession period
times = [0.0, jet.ne.modifier.period/3, 2*jet.ne.modifier.period/3]

images = map(times) do t
    cam = S.CameraZ(; xys=grid(SVector, x=range(2±extent, npx), y=range(0±extent, npx)), nz, ν, t)
    S.render(cam, jet, S.Intensity())
end


# ===== Create figure with 3 horizontal panels =====
fig = Figure(; size=(1200, 400))

# Compute common color limits
all_vals = vcat(vec.(images)...)
colorrange = maximum(all_vals) .* (0, 1)
colorscale = SymLog(1e-3*maximum(all_vals))

for (i, (img, t)) in enumerate(zip(images, times))
    ax = Axis(fig[1, i];
        title = "t = $(round(t, digits=1))",
        xlabel = i == 2 ? "x" : "",
        ylabel = i == 1 ? "y" : "",
        aspect = DataAspect()
    )
    
    heatmap!(ax, img; colorscale, colorrange)
end

Colorbar(fig[1, 4]; label="Intensity", scale=colorscale, colorrange)

save("rotating_beam_times.png", fig)
display(fig)

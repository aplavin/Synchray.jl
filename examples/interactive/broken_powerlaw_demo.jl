using DataManipulation
using MakieExtra; import GLMakie
using Unitful
using AxisKeysExtra, RectiGrids
using StructArrays
using AccessibleModels; using AccessibleModels: P
using IntervarrrSets
using DistributionsExtra
using Rotations
using OhMyThreads: tmap
using Metal
using SQLCollections, SQLite

MakieExtra.show_gl_icon_in_dock()

using Synchray
import Synchray as S


const ν₁ = 1.0
const ν₂ = 5.0

_tmap(f, X::KeyedArray) = @set AxisKeys.keyless_unname(X) = tmap(f, X)


sstate = let
    db = SQLite.DB(joinpath(@__DIR__, "_interactive_states.db"))
    SQLDictionary{String, Float64}(db, :broken_powerlaw_demo_state)
end


fig = Figure(size=(1600, 900))

params, = SliderGrid(fig[1:2,1][1,1], AccessibleModel((;
    img=(
        npix=P(discreterange(log, 16..500, length=50)),
        nz=P(discreterange(log, 2..300, length=50)),
        dynrange=P(LogUniform(1..1e5), 1e3),
    ),
    geom=(
        z=1..500,
        viewing_ang=P(LogUniform(0.1..80))u"°",
        opening_ang=P(LogUniform(0.1..80))u"°",
    ),
    phys=(
        γ=P(1.01..50, 5),
        B0=P(LogUniform(1e-10..1e-4), 1e-7),
        ne0=P(LogUniform(1e2..1e18), 1e6),
    ),
    electrons=(
        p_low=P(2:0.25:4, 2.5),
        p_high=P(2:0.25:4, 2.5),
        γ_break_inner=P(LogUniform(1..1000)),
        γ_break_outer=P(LogUniform(1..1000)),
    ),
), AccessibleModels.Auto()); state=sstate, rowgap=2)


# --- Two cameras at ν₁ and ν₂ ---

camera₁ = @lift S.CameraZ(;
    xys=grid(SVector,
            x=range(-5..30, length=2*$params.img.npix) |> LinRange,
            y=range(-5..5, length=$params.img.npix) |> LinRange),
    ν=ν₁, nz=$params.img.nz, t=0.0,
    mapfunc=_tmap)

camera₂ = @lift S.CameraZ(;
    xys=grid(SVector,
            x=range(-5..30, length=2*$params.img.npix) |> LinRange,
            y=range(-5..5, length=$params.img.npix) |> LinRange),
    ν=ν₂, nz=$params.img.nz, t=0.0,
    mapfunc=_tmap)


# --- Jet model ---

jet = @lift let
    axis = RotY($params.geom.viewing_ang) * SVector(0,0,1)

    region = S.EmissionRegion(;
        geometry = S.Geometries.Conical(; axis, φj=$params.geom.opening_ang, z=$params.geom.z),
        velocity = S.VelocitySpec(S.Directions.Radial(), S.gamma, S.Profiles.Transverse(Returns($params.phys.γ))),
        emission = S.SynchrotronEmission(;
            ne = S.Profiles.Axial(S.PowerLaw(-2; val0=$params.phys.ne0, s0=1.0)),
            B = S.BFieldSpec(
                S.Profiles.Axial(S.PowerLaw(-1; val0=$params.phys.B0, s0=1.0)),
                S.Directions.Scalar(),
                b -> S.FullyTangled(b),
            ),
            electrons = let
                γ_bi = $params.electrons.γ_break_inner
                γ_bo = $params.electrons.γ_break_outer
                S.BrokenPowerLaw(S.IsotropicPowerLawElectrons;
                    p_low=$params.electrons.p_low,
                    p_high=$params.electrons.p_high,
                    γ_break=S.Profiles.Transverse(η -> γ_bi + (γ_bo - γ_bi) * clamp(2η, 0, 1)),
                    Cj=1.0, Ca=1.0)
            end,
        ),
    )

    S.prepare_for_computations(region)
end::Any


# --- GPU preparation ---

camera₁_gpu = @lift @p $camera₁ S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))
camera₂_gpu = @lift @p $camera₂ S.to_float_type(Float32) @set(__.mapfunc = map) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))
jet_gpu = @lift S.to_float_type(Float32, $jet)::Any


# --- Render all four images ---

img_I₁ = @lift let
    (value, time) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera₁_gpu, $jet_gpu, S.Intensity())))
    @info "Intensity ν₁" time
    value
end

img_SI₁ = @lift let
    (value, time) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera₁_gpu, $jet_gpu, S.SpectralIndex())))
    @info "SpectralIndex ν₁" time
    value
end

img_I₂ = @lift let
    (value, time) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera₂_gpu, $jet_gpu, S.Intensity())))
    @info "Intensity ν₂" time
    value
end

img_SI₂ = @lift let
    (value, time) = @timed @modify(Array, AxisKeys.keyless_unname(S.render($camera₂_gpu, $jet_gpu, S.SpectralIndex())))
    @info "SpectralIndex ν₂" time
    value
end


# --- Interactive spectrum probe points ---

spec_points = Observable([(x=0., y=0.)])
ip = InteractivePoints(spec_points)
COLORS = Makie.wong_colors()

SPEC_FREQS = 10 .^ range(log10(0.1*ν₁), log10(10*ν₂), length=500)


# --- Display: intensity images ---

# Row 1: ν₁
ax_I₁, hm = image(fig[1,2][1,1], img_I₁;
    colormap=:turbo,
    colorscale=(@lift SymLog(1/$params.img.dynrange * maximum($img_I₁))),
    colorrange=(@lift (0, 1) .* maximum($img_I₁)),
    interpolate=false,
    axis=(; title="Intensity  ν=$ν₁"),
)
Colorbar(fig[1,2][1,2], hm, tickformat=EngTicks())

# Add interactive points on ν₁ image
axplot(scatter!, widgets=[ip])(ax_I₁,
    (@lift FPlot($spec_points, (@o _.x), (@o _.y)));
    color=(@lift [COLORS[mod1(i, length(COLORS))] for i in 1:length($spec_points)]),
    markersize=12, strokewidth=2, strokecolor=:white)

# Row 2: ν₂
ax_I₂, hm = image(fig[2,2][1,1], img_I₂;
    colormap=:turbo,
    colorscale=(@lift SymLog(1/$params.img.dynrange * maximum($img_I₂))),
    colorrange=(@lift (0, 1) .* maximum($img_I₂)),
    interpolate=false,
    axis=(; title="Intensity  ν=$ν₂"),
)
Colorbar(fig[2,2][1,2], hm, tickformat=EngTicks())

# Mirror points on ν₂ image (non-interactive)
scatter!(ax_I₂,
    (@lift [pt.x for pt in $spec_points]),
    (@lift [pt.y for pt in $spec_points]);
    color=(@lift [COLORS[mod1(i, length(COLORS))] for i in 1:length($spec_points)]),
    markersize=12, strokewidth=2, strokecolor=:white)

# ax_SI₁, hm = image(fig[1,3][1,1], img_SI₁;
#     colormap=:turbo, colorrange=(-2, 1),
#     interpolate=false,
#     axis=(; title="Spectral index  ν=$ν₁"),
# )
# Colorbar(fig[1,3][1,2], hm)

# ax_SI₂, hm = image(fig[2,3][1,1], img_SI₂;
#     colormap=:turbo, colorrange=(-2, 1),
#     interpolate=false,
#     axis=(; title="Spectral index  ν=$ν₂"),
# )
# Colorbar(fig[2,3][1,2], hm)

linkaxes!(ax_I₁, ax_I₂)
# linkaxes!(ax_I₁, ax_SI₁, ax_I₂, ax_SI₂)


# --- Spectrum computation ---

spectra = @lift let
    pts = $spec_points
    j = $jet
    nz = $camera₁.nz
    t = $camera₁.t
    map(pts) do pt
        Iν = map(SPEC_FREQS) do ν
            ray = S.RayZ(; x0=S.FourPosition(t, pt.x, pt.y, 0.0), k=ν, nz)
            S.render(ray, j, S.Intensity())
        end
        (; ν=collect(SPEC_FREQS), I=Iν)
    end
end


# --- Spectrum panel ---

ax_spec = Axis(fig[1:2, 1][2,1];
    xscale=log10, yscale=log10,
    xlabel="ν", ylabel="I(ν)",
    title="Spectrum at probe points")

vlines!(ax_spec, [ν₁, ν₂]; color=:gray, linestyle=:dash, linewidth=1)

on(spectra) do specs
    # Remove old spectrum lines (keep the vlines which was added first)
    while length(ax_spec.scene.plots) > 1
        delete!(ax_spec, ax_spec.scene.plots[end])
    end
    for (i, sp) in enumerate(specs)
        c = COLORS[mod1(i, length(COLORS))]
        valid = sp.I .> 0
        any(valid) && lines!(ax_spec, sp.ν[valid], sp.I[valid]; color=c, linewidth=2)
    end
end
notify(spectra)

fig |> display

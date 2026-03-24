import Pkg
Pkg.activate(@__DIR__)

using Synchray
import Synchray as S
using MakieExtra; import CairoMakie
using RectiGrids

outdir = joinpath(@__DIR__, "images")
isdir(outdir) || mkpath(outdir)


# RGB rendering frequencies (code units, peaked at ν₀=1 for green)
const νs_RGB = [0.7, 1.0, 1.4]

function render_rgb(cam_base, obj)
    map(νs_RGB) do ν
        cam = @set cam_base.ν = ν
        S.render(cam, obj)
    end
end

function assemble_rgb(channels)
    # normalize each channel to [0,1] using a shared max across all channels
    mx = maximum(maximum, channels)
    mx == 0 && return fill(RGBf(0,0,0), size(channels[1]))
    map(channels[1], channels[2], channels[3]) do r, g, b
        RGBf(r/mx, g/mx, b/mx)
    end
end


function colored_ellipsoid_doppler()
    cam_base = S.CameraZ(;
        xys = grid(SVector, x=range(-2.0..2.0, 256), y=range(-2.0..2.0, 256)),
        nz = 256,
        ν = 1.0,  # placeholder, overridden per channel
        t = 0.0,
    )

    # Natively green ellipsoid: peaked emission at ν₀ = ν_G = 1.0
    # Convention: +z = toward observer, so βz>0 → blueshift, βz<0 → redshift
    make_obj(βz; α_val) = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0, 5),
                S.FourVelocity(SVector(0.0, 0.0, βz)),
            ),
            sizes = SVector(1.0, 1.0, 1.0),
        ),
        emission = S.PeakedEmission(
            S = S.Profiles.Constant(1.0),
            α = S.Profiles.Constant(α_val),
            ν₀ = 1.0,  # peak emission at ν₀ = 1.0 (green)
            σ = 0.3,
        ),
    )

    velocity_cases = [
        (label="away (β=−0.3)",    βz=-0.3),
        (label="stationary",        βz=0.0),
        (label="toward (β=+0.3)",   βz=+0.3),
    ]

    opacity_cases = [
        (label="opaque",      α_val=5.0),
        (label="transparent",  α_val=0.3),
    ]

    fig = Figure(size=(900, 650))
    for (row, opacity) in enumerate(opacity_cases), (col, vel) in enumerate(velocity_cases)
        obj = make_obj(vel.βz; opacity.α_val)
        channels = render_rgb(cam_base, obj)
        rgb = assemble_rgb(channels)

        title = row == 1 ? vel.label : ""
        ylabel = col == 1 ? opacity.label : ""
        ax = Axis(fig[row, col]; title, ylabel, aspect=DataAspect())
        image!(ax, rgb)
        hidespines!(ax)
        hidedecorations!(ax; label=false)
    end

    resize_to_layout!(fig)
    save(joinpath(outdir, "colored_ellipsoid_doppler.png"), fig)
    fig
end

function combined_three_ellipsoids()
    cam_base = S.CameraZ(;
        xys = grid(SVector, x=range(-4.0..4.0, 512), y=range(-4.0..4.0, 512)),
        nz = 256,
        ν = 1.0,
        t = 0.0,
    )

    # Back (low z): large stationary elliptic blob, green
    back = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0, 0),
                S.FourVelocity(SVector(0.0, 0.0, 0.0)),
            ),
            sizes = SVector(4, 1.3, 2),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    # Middle: r=2 sphere moving away (βz < 0 → redshift), center at xy=(1,0)
    mid = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 1, 0, 7),
                S.FourVelocity(SVector(0.0, 0.0, -0.2)),
            ),
            sizes = SVector(2.0, 2.0, 2.0),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    # Front (high z): r=2 sphere moving toward us (βz > 0 → blueshift), center at xy=(-1,0)
    front = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, -1, 0, 11),
                S.FourVelocity(SVector(0.0, 0.0, 0.2)),
            ),
            sizes = SVector(2.0, 2.0, 2.0),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    combined = S.CombinedMedium(back, mid, front)

    channels = render_rgb(cam_base, combined)
    rgb = assemble_rgb(channels)

    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1]; title="CombinedMedium: 3 ellipsoids", aspect=DataAspect())
    image!(ax, rgb)
    hidespines!(ax)
    hidedecorations!(ax)
    resize_to_layout!(fig)
    save(joinpath(outdir, "combined_three_ellipsoids.png"), fig)
    fig
end

function assemble_rgb(channels, mx; gamma=1)
    mx == 0 && return fill(RGBf(0,0,0), size(channels[1]))
    map(channels[1], channels[2], channels[3]) do r, g, b
        RGBf(clamp(r/mx, 0, 1)^(1/gamma), clamp(g/mx, 0, 1)^(1/gamma), clamp(b/mx, 0, 1)^(1/gamma))
    end
end

function animated_flying_ellipsoid(; nframes=120)
    cam_base = S.CameraZ(;
        xys = grid(SVector, x=range(-5.0..5.0, 256), y=range(-5.0..5.0, 256)),
        nz = 256,
        ν = 1.0,
        t = 0.0,
    )

    flyer_emission = S.PeakedEmission(S=S.Profiles.Constant(10.0), α=S.Profiles.Constant(0.5), ν₀=1.1, σ=0.15)
    flyer_sizes = SVector(2.0, 1.0, 1.0)
    make_flyer(x0, β) = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(x0, S.FourVelocity(β)),
            sizes = flyer_sizes,
        ),
        emission = flyer_emission,
    )

    # Background: large stationary green ellipsoid (low z = back)
    background = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0, 3),
                S.FourVelocity(SVector(0.0, 0.0, 0.0)),
            ),
            sizes = SVector(3.0, 2.0, 1.5),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.1, σ=0.15),
    )

    β = 0.6
    ang = deg2rad(40)

    # Center (y=0): pure transverse — redshifted by transverse Doppler
    flyer_center = make_flyer(S.FourPosition(0, 0, 0, 8), SVector(β, 0.0, 0.0))
    # Top (y>0): angled toward observer — blueshifted (δ≈1.38)
    flyer_approaching = make_flyer(S.FourPosition(0, 0, 2.5, 8), SVector(β * cos(ang), 0.0, β * sin(ang)))
    # Bottom (y<0): angled away from observer — strongly redshifted (δ≈0.38)
    flyer_receding = make_flyer(S.FourPosition(0, 0, -2.5, 8), SVector(β * cos(ang), 0.0, -β * sin(ang)))

    combined = S.CombinedMedium(background, flyer_center, flyer_approaching, flyer_receding)

    ts = range(-16.0, -1.0, length=nframes)

    # Pre-render all frames to find global normalization
    all_channels = map(ts) do t
        cam = @set cam_base.t = t
        render_rgb(cam, combined)
    end
    mx = maximum(ch -> maximum(maximum, ch), all_channels)

    # Record animation with gamma tone mapping to compress dynamic range
    γ = 3
    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1]; aspect=DataAspect())
    hidespines!(ax); hidedecorations!(ax)
    rgb0 = assemble_rgb(first(all_channels), mx; gamma=γ)
    img_obs = Observable(rgb0)
    image!(ax, img_obs)

    record(fig, joinpath(outdir, "animated_flying_ellipsoid.mp4"), eachindex(ts); framerate=30) do i
        img_obs[] = assemble_rgb(all_channels[i], mx; gamma=γ)
    end
    fig
end

function combined_three_ellipsoids_perspective()
    cam_base = S.CameraPerspective(;
        photon_direction = SVector(0.0, 0.0, 1.0),
        origin = SVector(0.0, 0.0, 20.0),
        xys = grid(SVector, x=range(-0.5..0.5, 512), y=range(-0.5..0.5, 512)),
        nz = 256,
        ν = 1.0,
        t = 0.0,
    )

    # Back (low z): large stationary elliptic blob, green
    back = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0, 0),
                S.FourVelocity(SVector(0.0, 0.0, 0.0)),
            ),
            sizes = SVector(4, 1.3, 2),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    # Middle: r=2 sphere moving away (βz < 0 → redshift), center at xy=(1,0)
    mid = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 1, 0, 7),
                S.FourVelocity(SVector(0.0, 0.0, -0.2)),
            ),
            sizes = SVector(2.0, 2.0, 2.0),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    # Front (high z): r=2 sphere moving toward us (βz > 0 → blueshift), center at xy=(-1,0)
    front = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, -1, 0, 11),
                S.FourVelocity(SVector(0.0, 0.0, 0.2)),
            ),
            sizes = SVector(2.0, 2.0, 2.0),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    combined = S.CombinedMedium(back, mid, front)

    channels = render_rgb(cam_base, combined)
    rgb = assemble_rgb(channels)

    fig = Figure(size=(600, 600))
    ax = Axis(fig[1, 1]; title="CameraPerspective: 3 ellipsoids", aspect=DataAspect())
    image!(ax, rgb)
    hidespines!(ax)
    hidedecorations!(ax)
    resize_to_layout!(fig)
    save(joinpath(outdir, "combined_three_ellipsoids_perspective.png"), fig)
    fig
end

function perspective_xy_movers(; nframes=60)
    cam_base = S.CameraPerspective(;
        photon_direction = SVector(0.0, 0.0, 1.0),
        origin = SVector(0.0, 0.0, 9.0),
        xys = grid(SVector, x=range(-1..1, 256), y=range(-1..1, 256)),
        nz = 256,
        ν = 1.0,
        t = 0.0,
    )

    # Stationary background ellipsoid at z=0
    bg = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0, 0),
                S.FourVelocity(SVector(0.0, 0.0, 0.0)),
            ),
            sizes = SVector(3.0, 3.0, 1.5),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.5), ν₀=1.0, σ=0.15),
    )

    # Sphere moving in +x direction
    mover_x = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0, 0.5, 5),
                S.FourVelocity(SVector(0.9, 0.0, 0.0)),
            ),
            sizes = SVector(1.5, 1.5, 1.5),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.8), ν₀=1.0, σ=0.15),
    )

    # Sphere moving in -y direction
    mover_y = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, -0.5, 0, 3),
                S.FourVelocity(SVector(0.0, -0.9, 0.0)),
            ),
            sizes = SVector(1.5, 1.5, 1.5),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.8), ν₀=1.0, σ=0.15),
    )

    # Sphere moving diagonally in xy plane
    mover_xy = S.EmissionRegion(
        geometry = S.Geometries.Ellipsoid(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, 0.5, -0.5, 1),
                S.FourVelocity(SVector(-0.6, 0.6, 0.0)),
            ),
            sizes = SVector(1.5, 1.5, 1.5),
        ),
        emission = S.PeakedEmission(S=S.Profiles.Constant(1.0), α=S.Profiles.Constant(0.8), ν₀=1.0, σ=0.15),
    )

    combined = S.CombinedMedium(bg, mover_x, mover_y, mover_xy)

    # Side-by-side animation: SlowLight vs FastLight
    ts = range(-20.0, 20.0, length=nframes)
    γ = 6

    all_ch = map([S.SlowLight(), S.FastLight()]) do light_mode
        cam_l = @set cam_base.light = light_mode
        map(ts) do t
            cam = @set cam_l.t = t
            render_rgb(cam, combined)
        end
    end
    mx = maximum(chs -> maximum(ch -> maximum(maximum, ch), chs), all_ch)

    fig = Figure(size=(1200, 600))
    ax1 = Axis(fig[1, 1]; title="SlowLight", aspect=DataAspect())
    ax2 = Axis(fig[1, 2]; title="FastLight", aspect=DataAspect())
    for ax in (ax1, ax2); hidespines!(ax); hidedecorations!(ax; label=false); end

    obs_sl = Observable(assemble_rgb(first(all_ch[1]), mx; gamma=γ))
    obs_fl = Observable(assemble_rgb(first(all_ch[2]), mx; gamma=γ))
    image!(ax1, obs_sl)
    image!(ax2, obs_fl)

    # Save PNG at t=0 frame
    i_zero = argmin(abs.(ts))
    obs_sl[] = assemble_rgb(all_ch[1][i_zero], mx; gamma=γ)
    obs_fl[] = assemble_rgb(all_ch[2][i_zero], mx; gamma=γ)
    save(joinpath(outdir, "perspective_xy_movers.png"), fig)

    record(fig, joinpath(outdir, "perspective_xy_movers.mp4"), eachindex(ts); framerate=7) do i
        obs_sl[] = assemble_rgb(all_ch[1][i], mx; gamma=γ)
        obs_fl[] = assemble_rgb(all_ch[2][i], mx; gamma=γ)
    end
    fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    colored_ellipsoid_doppler()
    combined_three_ellipsoids()
    combined_three_ellipsoids_perspective()
    perspective_xy_movers()
    animated_flying_ellipsoid(; nframes=20)
end

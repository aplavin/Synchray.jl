include("helpers.jl")

using Random


function make_ball(; center_xyz, z, β=SVector(0.0, 0.0, 0.0), radius=1.0, S_val=1.0, α_val=0.5, ν₀=1.0, σ=0.15)
    S.EmissionRegion(
        geometry = S.Geometries.Ball(;
            center = S.Geometries.InertialWorldline(
                S.FourPosition(0, center_xyz..., z),
                S.FourVelocity(β),
            ),
            size = radius,
        ),
        emission = S.PeakedEmission(;
            S = S.Profiles.Constant(S_val),
            α = S.Profiles.Constant(α_val),
            ν₀, σ,
        ),
    )
end

function combined_medium_benchmark(; render_field, log_path, suffix="")
    log_section!(log_path, "combined_medium_benchmark")

    # --- Three objects (shared by tuple and vector panels) ---
    obj_back = make_ball(center_xyz=(0.0, 0.5), z=0.0, β=SVector(0.1, 0.0, 0.0), radius=2.0)
    obj_mid  = make_ball(center_xyz=(1.0, -0.3), z=5.0, β=SVector(0.0, 0.0, -0.15), radius=1.2, ν₀=0.9)
    obj_front = make_ball(center_xyz=(-0.8, 0.0), z=9.0, β=SVector(-0.1, 0.05, 0.1), radius=1.0, ν₀=1.1)

    cm_tuple = S.CombinedMedium(obj_back, obj_mid, obj_front)
    cm_vec3  = S.CombinedMedium([obj_back, obj_mid, obj_front])

    # --- 100 objects spread across the field of view ---
    rng = Random.MersenneTwister(42)
    n_objs = 300
    z_spacing = 1.0  # ensure non-overlapping along z
    objs_100 = map(1:n_objs) do i
        xy = SVector(randn(rng) * 2.5, randn(rng) * 2.5)
        z = (i - 1) * z_spacing
        β = SVector(0.1 * randn(rng), 0.1 * randn(rng), 0.05 * randn(rng))
        r = 0.2 + 0.15 * rand(rng)
        ν₀ = 0.8 + 0.4 * rand(rng)
        make_ball(; center_xyz=Tuple(xy), z, β, radius=r, ν₀)
    end
    cm_vec100 = S.CombinedMedium(objs_100)

    # --- Render all three ---
    img_tuple = render_field(cm_tuple; extent=4, ν=1.0, what=S.Intensity())
    img_vec3  = render_field(cm_vec3;  extent=4, ν=1.0, what=S.Intensity())
    img_vec100 = render_field(cm_vec100; extent=6, ν=1.0, what=S.Intensity())

    # --- Figure ---
    fig = Figure(size=(1200, 400))
    for (col, (img, title)) in enumerate([
        (img_tuple, "3 objects (Tuple)"),
        (img_vec3,  "3 objects (Vector)"),
        (img_vec100, "100 objects (Vector)"),
    ])
        ax = Axis(fig[1, col]; title, aspect=DataAspect())
        heatmap!(ax, img; colormap=:magma, colorscale=SymLog(1e-4 * maximum(img)))
        hidespines!(ax)
        hidedecorations!(ax)
    end
    resize_to_layout!(fig)
    save(joinpath(outdir, "combined_medium$(suffix).png"), fig)
    fig
end


function main(; suffix="")
    logs = init_script(@__FILE__)
    combined_medium_benchmark(; render_field=make_render_cpu(logs.log_cpu), log_path=logs.log_cpu, suffix)
end

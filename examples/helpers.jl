import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()

outdir = joinpath(@__DIR__, "images")
isdir(outdir) || mkpath(outdir)

const HAS_METAL = Sys.isapple()
@static if HAS_METAL
    using Metal
end

using Synchray
import Synchray as S
using MakieExtra; import CairoMakie
using RectiGrids
using AxisKeysExtra
using PyFormattedStrings
using Accessors


function init_script(script_file)
    timedir = joinpath(@__DIR__, "render_times")
    isdir(timedir) || mkpath(timedir)

    base = splitext(basename(script_file))[1]
    log_cpu = joinpath(timedir, "$(base).txt")
    open(log_cpu, "w") do _ end

    log_metal = joinpath(timedir, "$(base)_metal.txt")
    HAS_METAL && open(log_metal, "w") do _ end

    (; log_cpu, log_metal)
end

log_section!(log_path, name) = open(io -> write(io, "$name()\n"), log_path, "a")

function make_render_cpu(log_path)
    function(obj; extent, nz=256, npx=256, ν, t=0.0, what=S.Intensity(), adaptive_supersampling=false)
        cam = S.CameraZ(; xys=grid(SVector, x=range(0±extent, npx), y=range(0±extent, npx)), nz, ν, t)
        open(log_path, "a") do io
            redirect_stdout(io) do
                S.render(cam, obj, what; adaptive_supersampling)
                label = f"{nameof(typeof(obj)):27s}, {nameof(typeof(what)):15s}, {npx}px × {nz}"
                @time label S.render(cam, obj, what; adaptive_supersampling)
            end
        end
    end
end

function make_render_metal(log_path)
    function(obj; extent, nz=512, npx=512, ν, t=0.0, what=S.Intensity(), adaptive_supersampling=false)
        cam = S.CameraZ(; xys=grid(SVector, x=range(0±extent, npx), y=range(0±extent, npx)), nz, ν=Float32(ν), t)
        obj32 = S.to_float_type(Float32, obj)
        cam_gpu = @p cam S.to_float_type(Float32) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))
        open(log_path, "a") do io
            redirect_stdout(io) do
                @modify(Array, AxisKeys.keyless_unname(S.render(cam_gpu, obj32, what)))
                label = f"{nameof(typeof(obj)):27s}, {nameof(typeof(what)):15s}, {npx}px × {nz}"
                @time label @modify(Array, AxisKeys.keyless_unname(S.render(cam_gpu, obj32, what)))
            end
        end
    end
end


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


# RGB rendering helpers
const νs_RGB = [0.7, 1.0, 1.4]

function render_rgb(cam_base, obj; log_path=nothing)
    do_render() = map(νs_RGB) do ν
        cam = @set cam_base.ν = ν
        S.render(cam, obj)
    end
    if log_path !== nothing
        do_render()
        npx = size(cam_base.xys, 1)
        nz = cam_base.nz
        label = f"{nameof(typeof(obj)):27s},     RGB 3×render, {npx}px × {nz}"
        open(log_path, "a") do io
            redirect_stdout(io) do
                @time label do_render()
            end
        end
    else
        do_render()
    end
end

function render_rgb_metal(cam_base, obj; log_path=nothing)
    obj32 = S.to_float_type(Float32, obj)
    cam32_base = @p cam_base S.to_float_type(Float32) @modify(MtlArray, AxisKeys.keyless_unname(__.xys))
    do_render() = map(νs_RGB) do ν
        cam_gpu = @set cam32_base.ν = Float32(ν)
        @modify(Array, AxisKeys.keyless_unname(S.render(cam_gpu, obj32)))
    end
    if log_path !== nothing
        do_render()
        npx = size(cam_base.xys, 1)
        nz = cam_base.nz
        label = f"{nameof(typeof(obj)):27s},     RGB 3×render, {npx}px × {nz}"
        open(log_path, "a") do io
            redirect_stdout(io) do
                @time label do_render()
            end
        end
    else
        do_render()
    end
end

function assemble_rgb(channels)
    mx = maximum(maximum, channels)
    mx == 0 && return fill(RGBf(0,0,0), size(channels[1]))
    map(channels[1], channels[2], channels[3]) do r, g, b
        RGBf(r/mx, g/mx, b/mx)
    end
end

function assemble_rgb(channels, mx; gamma=1)
    mx == 0 && return fill(RGBf(0,0,0), size(channels[1]))
    map(channels[1], channels[2], channels[3]) do r, g, b
        RGBf(clamp(r/mx, 0, 1)^(1/gamma), clamp(g/mx, 0, 1)^(1/gamma), clamp(b/mx, 0, 1)^(1/gamma))
    end
end

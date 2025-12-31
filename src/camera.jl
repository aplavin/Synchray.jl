@kwdef struct CameraZ{Tν,Tt}
	xys
    nz::Int
    ν::Tν
    t::Tt
    mapfunc = map
end

struct RayZ{TX<:FourPosition, TK<:FourFrequency}
    x0::TX
    k::TK
    nz::Int
end
RayZ(; x0, k, nz::Int) = _ray_z(x0, k, nz)
_ray_z(x0::FourPosition, k::FourFrequency, nz::Int) = RayZ(x0, k, nz)
_ray_z(x0::FourPosition, k::Number, nz::Int) = RayZ(x0, photon_k(k, SVector(0, 0, 1)), nz)


photon_frequency(k::FourFrequency) = k.t
photon_frequency(ray::RayZ) = photon_frequency(ray.k)

Accessors.set(k::FourFrequency, ::typeof(photon_frequency), ν) = let
    ν0 = photon_frequency(k)
    return (ν / ν0) * k
end
Accessors.set(ray::RayZ, ::typeof(photon_frequency), ν) = @set photon_frequency(ray.k) = ν


struct Intensity end
struct OpticalDepth end
struct SpectralIndex end

render(ray::RayZ, obj::AbstractMedium, what=Intensity()) = integrate_ray(obj, ray, what)

_boundary_mask(img) = begin
    inside = map(x -> iszero(x) || isnan(x), img) |> Matrix
    IXs = CartesianIndices(inside)
    map(IXs) do I
        neighbors = filter(∈(IXs), I .+ CartesianIndices((-1:1, -1:1)))
        !allequal(inside[neighbors])
    end
end

_mean_samples(samples) = mean(skip(isnan, samples))
_mean_samples(samples::AbstractArray{<:Tuple}) = ntuple(i -> mean(s -> s[i], samples), length(first(samples)))

@unstable render(cam::CameraZ, obj::AbstractMedium, what=Intensity(); adaptive_supersampling=false) = let
    x0_base = FourPosition(cam.t, 0, 0, 0)
    k = photon_k(cam.ν, SVector(0, 0, 1))
    ray_base = RayZ(x0_base, k, cam.nz)
	img = cam.mapfunc(cam.xys) do xy
        ray = @set ray_base.x0 += FourPosition(0, xy..., 0)
		render(ray, obj, what)
	end

    adaptive_supersampling === false && return img
    n = adaptive_supersampling::Int
    @assert n ≥ 2

    steps = map(step, axiskeys(cam.xys))

    boundary = _boundary_mask(img)

    offs(n, d) = range(0 ± 0.7 * d, length=n)
    oxs, oys = offs.(n, steps)

    for I in findall(boundary)
        xy0 = cam.xys[I]
        samples = map(grid(SVector, oxs, oys)) do oxy
            xy = xy0 + oxy
            ray = @set ray_base.x0 += FourPosition(0, xy..., 0)
            render(ray, obj, what)
        end
        img[I] = _mean_samples(samples)
    end

    img
end

"""
    event_on_camera_ray(cam::CameraZ, r; t_obs=cam.t) -> FourPosition

Return the spacetime event `x4` on the camera's (orthographic, +z) null ray that passes
through spatial point `r = (x, y, z)`, at observer time `t_obs`.

This matches the internal `RayZ` convention used in transfer:

- The camera screen is at `z=0`.
- Rays propagate along `+z` with `x(z) = x0 + z*(1,0,0,1)`.

Therefore for a fixed observer time `t_obs` on the screen, the corresponding emission
event at depth `z` is `t = t_obs + z`.
"""
@inline event_on_camera_ray(cam::CameraZ, r::SVector{3}; t_obs=cam.t) =
    FourPosition(t_obs + r.z, r.x, r.y, r.z)

"""
    camera_ray_anchor(x4::FourPosition) -> FourPosition

Convert a lab-frame event `x4 = (t, x, y, z)` to the corresponding camera-ray anchor
event on the screen plane `z=0` for the (orthographic, +z) camera convention.

The returned event has:

- the same image-plane coordinates `(x, y)` (i.e. the pixel that sees `x4`),
- `z = 0`,
- the *camera time* (observer time on the screen) `t_cam = t - z`.
"""
@inline camera_ray_anchor(x4::FourPosition) = FourPosition(x4.t - x4.z, x4.x, x4.y, 0)

# XXX: ideally, should just be the below (split rays() vs render()),
# but somehow it results in a lot of allocations for map(mapview(...))

# @unstable rays(cam::CameraZ) = let
#     x0_base = FourPosition(cam.t, 0, 0, 0)
#     k = photon_k(cam.ν, SVector(0, 0, 1))
#     mapview(cam.xys) do xy
#         x0 = x0_base + FourPosition(0, xy..., 0)
#         RayZ(x0, k, cam.nz)
#     end
# end

# @unstable render(cam::CameraZ, obj::AbstractMedium, what=Intensity()) = begin
# 	cam.mapfunc(rays(cam)) do ray
# 		render(ray, obj, what)
# 	end
# end

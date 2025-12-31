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

@unstable render(cam::CameraZ, obj::AbstractMedium, what=Intensity()) = let
    x0_base = FourPosition(cam.t, 0, 0, 0)
    k = photon_k(cam.ν, SVector(0, 0, 1))
	cam.mapfunc(cam.xys) do xy
        x0 = x0_base + FourPosition(0, xy..., 0)
        ray = RayZ(x0, k, cam.nz)
		render(ray, obj, what)
	end
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

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
RayZ(; x0, k, nz::Int) = _rayz(x0, k, nz)
_rayz(x0::FourPosition, k::FourFrequency, nz::Int) = RayZ(x0, k, nz)
_rayz(x0::FourPosition, k::Number, nz::Int) = RayZ(x0, photon_k(k, SVector(0, 0, 1)), nz)

   
photon_frequency(k::FourFrequency) = k.t
photon_frequency(ray::RayZ) = photon_frequency(ray.k)

Accessors.set(k::FourFrequency, ::typeof(photon_frequency), ν) = let
    ν0 = photon_frequency(k)
    return (ν / ν0) * k
end
Accessors.set(ray::RayZ, ::typeof(photon_frequency), ν) = @set photon_frequency(ray.k) = ν


@unstable rays(cam::CameraZ) = let
    x0_base = FourPosition(cam.t, 0, 0, 0)
    k = photon_k(cam.ν, SVector(0, 0, 1))
    map(cam.xys) do xy
        x0 = x0_base + FourPosition(0, xy..., 0)
        RayZ(x0, k, cam.nz)
    end
end


struct Intensity end
struct OpticalDepth end
struct SpectralIndex end

render(ray::RayZ, obj::AbstractMedium, what=Intensity()) = integrate_ray(obj, ray, what)

@unstable render(cam::CameraZ, obj::AbstractMedium, what=Intensity()) = begin
	cam.mapfunc(rays(cam)) do ray
		render(ray, obj, what)
	end
end

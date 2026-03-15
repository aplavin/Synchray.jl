"""
	CameraOrtho(; look_direction, xys, nz, ν, t, origin=zeros(3), up=SVector(0,1,0), mapfunc=map)

Parallel-ray camera in flat spacetime.

All rays share the same propagation direction `n̂ = normalize(look_direction)` and the same
frequency `ν`. The camera defines a screen plane: an infinite 2D plane perpendicular to `n̂`
that passes through `origin`. Each pixel `(u, v)` in `xys` corresponds to a ray anchored at
`origin + u·e1 + v·e2`, propagating in the `n̂` direction.

# Fields
- `origin`: center of the screen plane — the spatial point where `(u, v) = (0, 0)` maps to.
  Only matters when the scene is not centered at the coordinate origin, or for `event_on_camera_ray`.
- `n`:  unit ray propagation direction (derived from `look_direction`).
- `e1`, `e2`: orthonormal screen-plane basis vectors (⊥ `n`, ⊥ each other). These also serve
  as the lab-frame polarization reference frame: Stokes Q is measured along `e1`. Derived from
  `up × n` (and then `n × e1`).
- `xys`: pixel coordinates — each element `(u, v)` is a 2D offset in the `(e1, e2)` basis.
- `nz`: number of integration steps per ray.
- `ν`: observing frequency.
- `t`: observer time — the time coordinate assigned to light arriving at the screen plane.
- `mapfunc`: mapping function applied over `xys` (default `map`; pass a GPU kernel for GPU rendering).
"""
struct CameraOrtho{Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf}
	origin::To
	n::Tn
	e1::Tn
	e2::Tn
	xys::Txys
    nz::Int
    ν::Tν
    t::Tt
    mapfunc::Tf
end

CameraOrtho(origin::To, n::Tn, e1::Tn, e2::Tn, xys::Txys, nz::Integer, ν::Tν, t::Tt, mapfunc::Tf=map) where {Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf} =
	CameraOrtho{Tν, Tt, To, Tn, Txys, Tf}(origin, n, e1, e2, xys, Int(nz), ν, t, mapfunc)

function CameraOrtho(; look_direction::SVector{3}, up=SVector(0, 1, 0), origin=zero(look_direction), xys, nz, ν, t, mapfunc=map)
	n = normalize(look_direction)
	e1 = normalize(cross(up, n))
	e2 = cross(n, e1)
	CameraOrtho(origin, n, e1, e2, xys, nz, ν, t, mapfunc)
end

"""Convenience constructor for the common Z-direction camera (screen at z=0, rays along +z)."""
function CameraZ(; xys, nz, ν, t, mapfunc=map)
	sample = first(xys)
	OT = float(eltype(sample))          # coordinate type (may have units)
	FT = typeof(float(one(OT)))         # dimensionless float type
	CameraOrtho(
		zero(SVector{3,OT}),
		SVector{3,FT}(0, 0, 1),
		SVector{3,FT}(1, 0, 0),
		SVector{3,FT}(0, 1, 0),
		xys, nz, ν, t, mapfunc
	)
end


struct Ray{TX<:FourPosition, TK<:FourFrequency, TE<:SVector{3}}
    x0::TX       # anchor 4-position on the ray
    k::TK        # 4-frequency: k = ν·(1, n̂)
    e1::TE       # polarization frame vector 1 (lab-frame spatial, ⊥ n̂)
    e2::TE       # polarization frame vector 2 (lab-frame spatial, ⊥ n̂, ⊥ e1)
    nz::Int
end

"""Unit spatial direction of ray propagation."""
direction3(ray::Ray) = direction3(ray.k)

"""Convenience constructor for Z-direction rays (n̂ = ẑ, screen basis x̂/ŷ)."""
function RayZ(x0::FourPosition, k::FourFrequency, nz::Int)
	T = float(eltype(k))
	Ray(x0, k, SVector{3,T}(1, 0, 0), SVector{3,T}(0, 1, 0), nz)
end
RayZ(; x0, k, nz::Int) = _ray_z(x0, k, nz)
_ray_z(x0::FourPosition, k::FourFrequency, nz::Int) = RayZ(x0, k, nz)
_ray_z(x0::FourPosition, k::Number, nz::Int) = RayZ(x0, photon_k(k, SVector(0, 0, 1)), nz)


@unstable ustrip(cam::CameraOrtho) = @p let
    cam
	@modify(_u_to_code(_, UCTX.L0), __.origin)
    @modify(_u_to_code(_, UCTX.L0), __.xys)
    @modify(_u_to_code(_, UCTX.ν0), __.ν)
    @modify(_u_to_code(_, UCTX.T0), __.t)
end


frequency(ray::Ray) = frequency(ray.k)
Accessors.set(ray::Ray, ::typeof(frequency), ν) = @set frequency(ray.k) = ν


struct Intensity end
struct IntensityIQU end
struct OpticalDepth end
struct SpectralIndex end

render(ray::Ray, obj::AbstractMedium, what=Intensity()) = integrate_ray(obj, ray, what)

@unstable _boundary_mask(img) = begin
    inside = map(x -> iszero(x) || isnan(x), img) |> Matrix
    IXs = CartesianIndices(inside)
    map(IXs) do I
        neighbors = filter(∈(IXs), I .+ CartesianIndices((-1:1, -1:1)))
        !allequal(inside[neighbors])
    end
end

_mean_samples(samples) = mean(skip(isnan, samples))
_mean_samples(samples::AbstractArray{<:Tuple}) = ntuple(i -> mean(s -> s[i], samples), length(first(samples)))

@unstable render(cam::CameraOrtho, obj::AbstractMedium, what=Intensity(); adaptive_supersampling=false) = let
    k = photon_k(cam.ν, cam.n)
    ray_base = Ray(FourPosition(cam.t, cam.origin), k, cam.e1, cam.e2, cam.nz)
	img = cam.mapfunc(cam.xys) do uv
        offset = uv[1] * cam.e1 + uv[2] * cam.e2
        ray = @set ray_base.x0 += FourPosition(0, offset)
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
        uv0 = cam.xys[I]
        samples = map(grid(SVector, oxs, oys)) do ouv
            uv = uv0 + ouv
            offset = uv[1] * cam.e1 + uv[2] * cam.e2
            ray = @set ray_base.x0 += FourPosition(0, offset)
            render(ray, obj, what)
        end
        img[I] = _mean_samples(samples)
    end

    img
end

"""
    event_on_camera_ray(cam::CameraOrtho, r; t_obs=cam.t) -> FourPosition

Return the spacetime event on the camera's null ray that passes through spatial point
`r = (x, y, z)`, at observer time `t_obs`.

The ray path is `x4(s) = x0 + s·(1, n̂)`, so `t = t_obs + s` where `s = dot(r - origin, n̂)`.
"""
@inline event_on_camera_ray(cam::CameraOrtho, r::SVector{3}; t_obs=cam.t) = let
    s = dot(r - cam.origin, cam.n)
    FourPosition(t_obs + s, r)
end

"""
    camera_ray_anchor(cam::CameraOrtho, x4) -> FourPosition

Convert a lab-frame event `x4` to the camera-ray anchor on the screen plane.

Returns the event with:
- spatial position projected onto the screen plane (perpendicular to n̂)
- camera time `t_cam = t - s` where `s = dot(r - origin, n̂)`
"""
@inline camera_ray_anchor(cam::CameraOrtho, x4::FourPosition) = let
    r = @swiz x4.xyz
    s = dot(r - cam.origin, cam.n)
    r_screen = r - s * cam.n
    FourPosition(x4.t - s, r_screen)
end

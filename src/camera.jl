"""Light propagation mode: controls the ray advancement direction through spacetime."""
abstract type AbstractLightMode end

"""Slow light: finite light speed (c = 1). Ray advances along the null direction `(1, n̂)`.
The time coordinate at depth s is `t = t_obs + s`. This is the default, physically correct mode."""
struct SlowLight <: AbstractLightMode end

"""Fast light: infinite light speed. Ray advances along the spatial direction `(0, n̂)`.
All events along the ray are at observer time `t = t_obs`."""
struct FastLight <: AbstractLightMode end

@inline _direction4(::SlowLight, n̂::SVector{3}) = FourPosition(one(eltype(n̂)), n̂)
@inline _direction4(::FastLight, n̂::SVector{3}) = FourPosition(zero(eltype(n̂)), n̂)


"""
	CameraOrtho(; photon_direction, xys, nz, ν, t, origin=zeros(3), up=SVector(0,1,0), light=SlowLight(), mapfunc=map)

Parallel-ray camera in flat spacetime.

All rays share the same propagation direction `n̂ = normalize(photon_direction)` and the same
frequency `ν`. The camera defines a screen plane: an infinite 2D plane perpendicular to `n̂`
that passes through `origin`. Each pixel `(u, v)` in `xys` corresponds to a ray anchored at
`origin + u·e1 + v·e2`, propagating in the `n̂` direction.

# Fields
- `origin`: center of the screen plane — the spatial point where `(u, v) = (0, 0)` maps to.
  Only matters when the scene is not centered at the coordinate origin, or for `event_on_camera_ray`.
- `n`:  unit photon propagation direction (derived from `photon_direction`).
- `e1`, `e2`: orthonormal screen-plane basis vectors (⊥ `n`, ⊥ each other). These also serve
  as the lab-frame polarization reference frame: Stokes Q is measured along `e1`. Derived from
  `up × n` (and then `n × e1`).
- `xys`: pixel coordinates — each element `(u, v)` is a 2D offset in the `(e1, e2)` basis.
- `nz`: number of integration steps per ray.
- `ν`: observing frequency.
- `t`: observer time — the time coordinate assigned to light arriving at the screen plane.
- `light`: light propagation mode (`SlowLight()` or `FastLight()`).
- `mapfunc`: mapping function applied over `xys` (default `map`; pass a GPU kernel for GPU rendering).
"""
struct CameraOrtho{Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf, TL<:AbstractLightMode}
	origin::To
	n::Tn
	e1::Tn
	e2::Tn
	xys::Txys
    nz::Int
    ν::Tν
    t::Tt
    light::TL
    mapfunc::Tf
end

CameraOrtho(origin::To, n::Tn, e1::Tn, e2::Tn, xys::Txys, nz::Integer, ν::Tν, t::Tt, light::TL=SlowLight(), mapfunc::Tf=map) where {Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf, TL<:AbstractLightMode} =
	CameraOrtho{Tν, Tt, To, Tn, Txys, Tf, TL}(origin, n, e1, e2, xys, Int(nz), ν, t, light, mapfunc)

function CameraOrtho(; photon_direction::SVector{3}, up=SVector(0, 1, 0), origin=zero(photon_direction), xys, nz, ν, t, light=SlowLight(), mapfunc=map)
	n = normalize(photon_direction)
	e1 = normalize(cross(up, n))
	e2 = cross(n, e1)
	CameraOrtho(origin, n, e1, e2, xys, nz, ν, t, light, mapfunc)
end

"""Convenience constructor for the common Z-direction camera (screen at z=0, photon direction +z)."""
@unstable function CameraZ(; xys, nz, ν, t, light=SlowLight(), mapfunc=map)
	sample = first(xys)
	OT = float(eltype(sample))          # coordinate type (may have units)
	FT = typeof(float(one(OT)))         # dimensionless float type
	CameraOrtho(
		zero(SVector{3,OT}),
		SVector{3,FT}(0, 0, 1),
		SVector{3,FT}(1, 0, 0),
		SVector{3,FT}(0, 1, 0),
		xys, nz, ν, t, light, mapfunc
	)
end


"""
	CameraPerspective(; photon_direction, origin, xys, nz, ν, t, up=SVector(0,1,0), light=SlowLight(), mapfunc=map)

Perspective (pinhole) camera in flat spacetime.

All rays originate from `origin` and fan out. The on-axis photon propagates in direction
`n̂ = normalize(photon_direction)`. Each pixel `(u, v)` in `xys` tilts the photon direction
away from `n̂`: the photon direction is `normalize(n̂ − u·e1 − v·e2)` (the sign inversion
is the standard pinhole/perspective effect). The `xys` values are tangents of the angle
from the optical axis, so `xys ∈ [-1, 1]` gives a 90° field of view.

# Fields
Same as `CameraOrtho` — `origin`, `n`, `e1`, `e2`, `xys`, `nz`, `ν`, `t`, `light`, `mapfunc` —
but `origin` is the camera (pinhole) position and `xys` are angular-tangent pixel coordinates.
"""
struct CameraPerspective{Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf, TL<:AbstractLightMode}
	origin::To
	n::Tn
	e1::Tn
	e2::Tn
	xys::Txys
	nz::Int
	ν::Tν
	t::Tt
	light::TL
	mapfunc::Tf
end

CameraPerspective(origin::To, n::Tn, e1::Tn, e2::Tn, xys::Txys, nz::Integer, ν::Tν, t::Tt, light::TL=SlowLight(), mapfunc::Tf=map) where {Tν, Tt, To<:SVector{3}, Tn<:SVector{3}, Txys, Tf, TL<:AbstractLightMode} =
	CameraPerspective{Tν, Tt, To, Tn, Txys, Tf, TL}(origin, n, e1, e2, xys, Int(nz), ν, t, light, mapfunc)

function CameraPerspective(; photon_direction::SVector{3}, up=SVector(0, 1, 0), origin=zero(photon_direction), xys, nz, ν, t, light=SlowLight(), mapfunc=map)
	n = normalize(photon_direction)
	e1 = normalize(cross(up, n))
	e2 = cross(n, e1)
	CameraPerspective(origin, n, e1, e2, xys, nz, ν, t, light, mapfunc)
end


struct Ray{TX<:FourPosition, TK<:FourFrequency, TE<:SVector{3}, TL<:AbstractLightMode, TS}
    x0::TX       # anchor 4-position on the ray
    k::TK        # 4-frequency: k = ν·(1, n̂)
    e1::TE       # polarization frame vector 1 (lab-frame spatial, ⊥ n̂)
    e2::TE       # polarization frame vector 2 (lab-frame spatial, ⊥ n̂, ⊥ e1)
    nz::Int
    light::TL
    s_max::TS    # maximum valid ray parameter (nothing = unbounded, 0 = half-ray for perspective)
end

"""Construct a `Ray` defaulting to `SlowLight()` and unbounded `s_max`."""
Ray(x0::FourPosition, k::FourFrequency, e1::SVector{3}, e2::SVector{3}, nz::Int) =
    Ray(x0, k, e1, e2, nz, SlowLight(), nothing)

"""Construct a `Ray` with given light mode and unbounded `s_max`."""
Ray(x0::FourPosition, k::FourFrequency, e1::SVector{3}, e2::SVector{3}, nz::Int, light::AbstractLightMode) =
    Ray(x0, k, e1, e2, nz, light, nothing)

"""Unit spatial direction of ray propagation."""
@inline direction3(ray::Ray) = direction3(ray.k)

"""
    direction4(ray::Ray) -> FourPosition

Ray path 4-direction: the spacetime advancement direction `d` such that `x(s) = x₀ + s·d`.

For `SlowLight`: `d = (1, n̂)` — null geodesic, same as `direction4(ray.k)`.
For `FastLight`: `d = (0, n̂)` — purely spatial, all events at observer time.
"""
@inline direction4(ray::Ray) = _direction4(ray.light, direction3(ray))

"""Convenience constructor for Z-direction rays (n̂ = ẑ, screen basis x̂/ŷ)."""
function RayZ(x0::FourPosition, k::FourFrequency, nz::Int; light=SlowLight())
	T = float(eltype(k))
	Ray(x0, k, SVector{3,T}(1, 0, 0), SVector{3,T}(0, 1, 0), nz, light)
end
RayZ(; x0, k, nz::Int, light=SlowLight()) = _ray_z(x0, k, nz, light)
_ray_z(x0::FourPosition, k::FourFrequency, nz::Int, light=SlowLight()) = RayZ(x0, k, nz; light)
_ray_z(x0::FourPosition, k::Number, nz::Int, light=SlowLight()) = RayZ(x0, photon_k(k, SVector(0, 0, 1)), nz; light)

@unstable ustrip(cam::CameraOrtho) = @p let
    cam
	@modify(_u_to_code(_, UCTX.L0), __.origin)
    @modify(_u_to_code(_, UCTX.L0), __.xys)
    @modify(_u_to_code(_, UCTX.ν0), __.ν)
    @modify(_u_to_code(_, UCTX.T0), __.t)
end

@unstable ustrip(cam::CameraPerspective) = @p let
    cam
	@modify(_u_to_code(_, UCTX.L0), __.origin)
    @modify(_u_to_code(_, UCTX.ν0), __.ν)
    @modify(_u_to_code(_, UCTX.T0), __.t)
end


"""Compute ray-geometry intersection clipped to the ray's valid parameter domain (`s ≤ ray.s_max`)."""
@inline z_interval_clipped(obj, ray::Ray) = _clip_z_interval(z_interval(obj, ray), ray.s_max)
@inline _clip_z_interval(seg, ::Nothing) = seg
@inline _clip_z_interval(seg, s_max) = leftendpoint(seg) .. min(rightendpoint(seg), s_max)

frequency(ray::Ray) = frequency(ray.k)
Accessors.set(ray::Ray, ::typeof(frequency), ν) = @set frequency(ray.k) = ν


struct Intensity end
struct IntensityIQU end
struct OpticalDepth end
struct SpectralIndex end

render(ray::Ray, obj::AbstractMedium, what=Intensity()) = integrate_ray(obj, ray, what)

@inline _orthographic_ray_base(cam::CameraOrtho) =
	Ray(FourPosition(cam.t, cam.origin), photon_k(cam.ν, cam.n), cam.e1, cam.e2, cam.nz, cam.light)

@inline _camera_ray(cam::CameraOrtho, uv, ray_base::Ray) = begin
	offset = uv[1] * cam.e1 + uv[2] * cam.e2
	@set ray_base.x0 += FourPosition(0, offset)
end

"""
	camera_ray(cam, uv) -> Ray

Construct the single camera ray for pixel coordinates `uv`.
"""
@inline camera_ray(cam::CameraOrtho, uv) = _camera_ray(cam, uv, _orthographic_ray_base(cam))

@inline camera_ray(cam::CameraPerspective, uv) = begin
	(; e1, e2, n) = cam
	# n is the photon direction. Pixel offsets tilt the photon away from n (pinhole inversion).
	photon_unnorm = n - uv[1] * e1 - uv[2] * e2
	photon_dir = photon_unnorm / norm(photon_unnorm)
	# Per-ray polarization basis: project camera e1 perpendicular to photon direction
	e1_proj = e1 - dot(e1, photon_dir) * photon_dir
	e1_ray = e1_proj / norm(e1_proj)
	e2_ray = cross(photon_dir, e1_ray)
	k = photon_k(cam.ν, photon_dir)
	# s_max = 0: scene is at negative s, behind-camera is at positive s
	Ray(FourPosition(cam.t, cam.origin), k, e1_ray, e2_ray, cam.nz, cam.light, zero(float(eltype(cam.origin))))
end

@unstable _boundary_mask(img) = begin
    inside = map(x -> iszero(x) || isnan(x), img) |> Matrix
    IXs = CartesianIndices(inside)
    map(IXs) do I
        @p let
            CartesianIndices((-1:1, -1:1))
            Iterators.map(I + _)
            Iterators.filter(∈(IXs))
            any(inside[_] != inside[I])
        end
    end
end

_mean_samples(samples) = mean(skip(isnan, samples))
_mean_samples(samples::AbstractArray{<:Tuple}) = ntuple(i -> mean(s -> s[i], samples), length(first(samples)))

@unstable render(cam::CameraOrtho, obj::AbstractMedium, what=Intensity(); adaptive_supersampling=false) = let
	ray_base = _orthographic_ray_base(cam)
	# strip xys so the closure captures only isbits fields (xys may be an array)
	cam_noxys = @set cam.xys = nothing
	img = cam.mapfunc(cam.xys) do uv
		ray = _camera_ray(cam_noxys, uv, ray_base)
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
            ray = _camera_ray(cam_noxys, uv, ray_base)
            render(ray, obj, what)
        end
        img[I] = _mean_samples(samples)
    end

    img
end

@unstable render(cam::CameraPerspective, obj::AbstractMedium, what=Intensity(); adaptive_supersampling=false) = let
	# strip xys so the closure captures only isbits fields (xys may be an array)
	cam_noxys = @set cam.xys = nothing
	img = cam.mapfunc(cam.xys) do uv
		render(camera_ray(cam_noxys, uv), obj, what)
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
            render(camera_ray(cam, uv), obj, what)
        end
        img[I] = _mean_samples(samples)
    end

    img
end

"""
    event_on_camera_ray(cam::CameraOrtho, r; t_obs=cam.t) -> FourPosition

Return the spacetime event on the camera ray that passes through spatial point
`r = (x, y, z)`, at observer time `t_obs`.

For `SlowLight`: `t = t_obs + s` where `s = dot(r - origin, n̂)` (null ray).
For `FastLight`: `t = t_obs` (all events simultaneous).
"""
@inline event_on_camera_ray(cam::CameraOrtho, r::SVector{3}; t_obs=cam.t) =
    _event_on_camera_ray(cam.light, cam, r, t_obs)

@inline _event_on_camera_ray(::SlowLight, cam, r, t_obs) = let
    s = dot(r - cam.origin, cam.n)
    FourPosition(t_obs + s, r)
end
@inline _event_on_camera_ray(::FastLight, cam, r, t_obs) = FourPosition(t_obs, r)

"""
    camera_ray_anchor(cam::CameraOrtho, x4) -> FourPosition

Convert a lab-frame event `x4` to the camera-ray anchor on the screen plane.

Returns the event with:
- spatial position projected onto the screen plane (perpendicular to n̂)
- camera time: `t_cam = t - s` for `SlowLight`, `t_cam = t` for `FastLight`
"""
@inline camera_ray_anchor(cam::CameraOrtho, x4::FourPosition) =
    _camera_ray_anchor(cam.light, cam, x4)

@inline _camera_ray_anchor(::SlowLight, cam, x4) = let
    r = @swiz x4.xyz
    s = dot(r - cam.origin, cam.n)
    r_screen = r - s * cam.n
    FourPosition(x4.t - s, r_screen)
end
@inline _camera_ray_anchor(::FastLight, cam, x4) = let
    r = @swiz x4.xyz
    s = dot(r - cam.origin, cam.n)
    r_screen = r - s * cam.n
    FourPosition(x4.t, r_screen)
end

"""
    image_coordinates(cam, r::SVector{3}) -> SVector{2}
    image_coordinates(cam, x4::FourPosition) -> SVector{2}

Project a spatial position (or spacetime event) onto the camera image plane,
returning `(u, v)` pixel coordinates in the `(e1, e2)` basis.

For `CameraOrtho`: orthographic projection (drop depth component).
For `CameraPerspective`: pinhole projection (divide by depth along optical axis).
"""
@inline function image_coordinates(cam::CameraOrtho, r::SVector{3})
    dr = r - cam.origin
    SA[dot(dr, cam.e1), dot(dr, cam.e2)]
end

@inline function image_coordinates(cam::CameraPerspective, r::SVector{3})
    dir = normalize(r - cam.origin)
    SA[-dot(dir, cam.e1) / dot(dir, cam.n), -dot(dir, cam.e2) / dot(dir, cam.n)]
end

@inline image_coordinates(cam, x4::FourPosition) = image_coordinates(cam, SVector(@swiz x4.xyz))

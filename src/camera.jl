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
	CameraOrtho(; look_direction, xys, nz, ν, t, origin=zeros(3), up=SVector(0,1,0), light=SlowLight(), mapfunc=map)

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

function CameraOrtho(; look_direction::SVector{3}, up=SVector(0, 1, 0), origin=zero(look_direction), xys, nz, ν, t, light=SlowLight(), mapfunc=map)
	n = normalize(look_direction)
	e1 = normalize(cross(up, n))
	e2 = cross(n, e1)
	CameraOrtho(origin, n, e1, e2, xys, nz, ν, t, light, mapfunc)
end

"""Convenience constructor for the common Z-direction camera (screen at z=0, rays along +z)."""
function CameraZ(; xys, nz, ν, t, light=SlowLight(), mapfunc=map)
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
	CameraPerspective(; look_direction, origin, xys, nz, ν, t, up=SVector(0,1,0), light=SlowLight(), mapfunc=map)

Perspective (pinhole) camera in flat spacetime.

All rays originate from `origin` and fan out, each pixel `(u, v)` in `xys` determining a
unique ray direction `normalize(n̂ + u·e1 + v·e2)`.  The `xys` values are tangents of the
angle from the optical axis, so `xys ∈ [-1, 1]` gives a 90° field of view.

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

function CameraPerspective(; look_direction::SVector{3}, up=SVector(0, 1, 0), origin=zero(look_direction), xys, nz, ν, t, light=SlowLight(), mapfunc=map)
	n = normalize(look_direction)
	e1 = normalize(cross(up, n))
	e2 = cross(n, e1)
	CameraPerspective(origin, n, e1, e2, xys, nz, ν, t, light, mapfunc)
end


struct Ray{TX<:FourPosition, TK<:FourFrequency, TE<:SVector{3}, TL<:AbstractLightMode}
    x0::TX       # anchor 4-position on the ray
    k::TK        # 4-frequency: k = ν·(1, n̂)
    e1::TE       # polarization frame vector 1 (lab-frame spatial, ⊥ n̂)
    e2::TE       # polarization frame vector 2 (lab-frame spatial, ⊥ n̂, ⊥ e1)
    nz::Int
    light::TL
end

"""Construct a `Ray` defaulting to `SlowLight()`."""
Ray(x0::FourPosition, k::FourFrequency, e1::SVector{3}, e2::SVector{3}, nz::Int) =
    Ray(x0, k, e1, e2, nz, SlowLight())

"""Unit spatial direction of ray propagation."""
direction3(ray::Ray) = direction3(ray.k)

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

"""Compute a unit vector perpendicular to `n` (preferring to stay close to `hint`)."""
@inline _perpendicular_basis_vector(n::SVector{3}, hint::SVector{3}) = begin
	v = cross(hint, n)
	nv = norm(v)
	if nv < 1e-8  # hint parallel to n; pick arbitrary perpendicular
		v = cross(SVector(0, 0, 1), n)
		nv = norm(v)
		if nv < 1e-8
			v = cross(SVector(1, 0, 0), n)
			nv = norm(v)
		end
	end
	v / nv
end


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
    (; e1, e2) = cam
    k = photon_k(cam.ν, cam.n)
    ray_base = Ray(FourPosition(cam.t, cam.origin), k, e1, e2, cam.nz, cam.light)
	img = cam.mapfunc(cam.xys) do uv
        offset = uv[1] * e1 + uv[2] * e2
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
            offset = uv[1] * e1 + uv[2] * e2
            ray = @set ray_base.x0 += FourPosition(0, offset)
            render(ray, obj, what)
        end
        img[I] = _mean_samples(samples)
    end

    img
end

@inline _perspective_ray(cam::CameraPerspective, uv) = begin
    (; e1, e2, n) = cam
    dir_unnorm = n + uv[1] * e1 + uv[2] * e2
    dir = dir_unnorm / norm(dir_unnorm)
    # Photon propagates from scene toward camera (opposite to pixel direction).
    # Using -dir for k gives correct Doppler shifts and SlowLight retardation.
    photon_dir = -dir
    # Per-ray polarization basis: project camera e1 perpendicular to photon direction
    e1_proj = e1 - dot(e1, photon_dir) * photon_dir
    e1_ray = e1_proj / norm(e1_proj)
    e2_ray = cross(photon_dir, e1_ray)
    k = photon_k(cam.ν, photon_dir)
    Ray(FourPosition(cam.t, cam.origin), k, e1_ray, e2_ray, cam.nz, cam.light)
end

@unstable render(cam::CameraPerspective, obj::AbstractMedium, what=Intensity(); adaptive_supersampling=false) = let
	img = cam.mapfunc(cam.xys) do uv
		render(_perspective_ray(cam, uv), obj, what)
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
            render(_perspective_ray(cam, uv), obj, what)
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


# ============================================================================
# GR Camera — gravitational lensing via precomputed deflection map
# ============================================================================

"""
	CameraGR(; camera, deflection, bh_position=zeros(3), bh_rg=1.0)

GR-lensed camera wrapping a flat-space `CameraOrtho`.

- `camera`: the underlying `CameraOrtho` (defines pixel grid, ν, t, nz, direction).
- `deflection`: array with the same shape as `camera.xys`, where each element is
  either a `Ray` (outgoing ray after BH deflection) or `nothing` (captured by BH).
- `bh_position`: spatial position of the black hole (default: origin).
- `bh_rg`: gravitational radius GM/c² (default: 1). Sets the coordinate scale for the
  deflection map — Krang computes in units of r_g, so outgoing ray coordinates are
  `bh_position + bh_rg * krang_coords`.

The deflection map is precomputed from the BH metric (see `compute_deflection_map`
in the Krang extension) and is reusable across emission models.
"""
@kwdef struct CameraGR{Tcam<:CameraOrtho, Tmap}
	camera::Tcam
	deflection::Tmap
	bh_position::SVector{3,Float64} = zero(SVector{3,Float64})
	bh_rg::Float64 = 1.0
end

@unstable render(cam::CameraGR, obj::AbstractMedium, what=Intensity()) = let
	(; camera, deflection, bh_position) = cam
	(; e1, e2) = camera
	k = photon_k(camera.ν, camera.n)
	ray_base = Ray(FourPosition(camera.t, camera.origin), k, e1, e2, camera.nz, camera.light)

	camera.mapfunc(camera.xys, deflection) do uv, ray_out
		offset = uv[1] * e1 + uv[2] * e2
		ray_in = @set ray_base.x0 += FourPosition(0, offset)

		_render_gr_pixel(obj, ray_in, ray_out, bh_position, what)
	end
end

"""s-parameter of closest approach to a point along a ray."""
@inline _s_closest_to_point(ray::Ray, point::SVector{3}) =
	dot(point - @swiz(ray.x0.xyz), direction3(ray))

"""Integrate through a medium (or CombinedMedium) with s-clipping, updating accumulator."""
@inline _integrate_through_clipped(acc, obj::AbstractMedium, ray, clip) = begin
	seg = z_interval(obj, ray) ∩ clip
	_integrate_segment(acc, obj, ray, seg)
end

@inline _integrate_through_clipped(acc, cm::CombinedMedium, ray, clip) =
	_integrate_clipped_recursive(acc, cm.objects, ray, clip)

@inline _integrate_clipped_recursive(acc, ::Tuple{}, ray, clip) = acc
@inline _integrate_clipped_recursive(acc, objs::Tuple, ray, clip) = begin
	obj = first(objs)
	seg = z_interval(obj, ray) ∩ clip
	acc = _integrate_segment(acc, obj, ray, seg)
	_integrate_clipped_recursive(acc, Base.tail(objs), ray, clip)
end

@inline _render_gr_pixel(obj, ray_in, ray_out::Nothing, bh_pos, what) = begin
	ν = frequency(ray_in)
	s_bh = _s_closest_to_point(ray_in, bh_pos)
	acc = _gr_init_acc(obj, ray_in, what)
	acc = _integrate_through_clipped(acc, obj, ray_in, typeof(s_bh)(-Inf) .. s_bh)
	_postprocess_acc(acc, ν, what)
end

@inline _render_gr_pixel(obj, ray_in, ray_out::Ray, bh_pos, what) = begin
	ν = frequency(ray_in)
	acc = _gr_init_acc(obj, ray_in, what)

	# Photon path: outgoing (far/background) → BH → incoming (near/foreground) → observer.
	# Integrate background first so foreground can attenuate it.

	# 1. Outgoing segment (background), clipped to outgoing half past BH
	s_bh_out = _s_closest_to_point(ray_out, bh_pos)
	acc = _integrate_through_clipped(acc, obj, ray_out, s_bh_out .. typeof(s_bh_out)(Inf))

	# 2. Incoming segment (foreground, attenuates background), clipped before BH
	s_bh_in = _s_closest_to_point(ray_in, bh_pos)
	acc = _integrate_through_clipped(acc, obj, ray_in, typeof(s_bh_in)(-Inf) .. s_bh_in)

	_postprocess_acc(acc, ν, what)
end

@inline _gr_init_acc(obj::AbstractMedium, ray, what) = begin
	seg = z_interval(obj, ray)
	k1 = direction4(ray)
	s0 = leftendpoint(seg)
	_integrate_ray_step(_init_acc(typeof(what), frequency(ray)), obj, ray.x0 + s0 * k1, ray, zero(float(s0)))
end

@inline _gr_init_acc(cm::CombinedMedium, ray, what) = _gr_init_acc(first(cm.objects), ray, what)

"""
Module for geometry types and coordinate transformations.

Geometries define the shape, boundaries, and natural coordinate systems
for emission regions.
"""
module Geometries

using ..Synchray
import ..Synchray as S

"""
    AbstractGeometry

Abstract type for emission region geometries.

Required interface:
- `z_interval(geom, ray)` - ray intersection interval
- `is_inside(geom, x4)` - containment test
- `natural_coords(geom, x4)` - natural coordinates as NamedTuple
- `geometry_axis(geom)` - primary axis direction
"""
abstract type AbstractGeometry end

# ============================================================================
# Conical Geometry
# ============================================================================

"""Helper struct for caching trig values"""
struct AngleTrigCached{T}
	sin::T
	cos::T
	tan::T
end

Base.sin(x::AngleTrigCached) = x.sin
Base.cos(x::AngleTrigCached) = x.cos
Base.tan(x::AngleTrigCached) = x.tan
Base.sincos(x::AngleTrigCached) = (x.sin, x.cos)
AngleTrigCached_fromangle(φ) = AngleTrigCached(sin(φ), cos(φ), tan(φ))

"""
    Conical{TR, Tφ, Tz}

Conical geometry with rotation matrix, half-opening angle, and axial extent.

# Fields
- `R_local_to_lab::TR`: 3×3 rotation matrix that maps local jet frame to lab frame (columns are local basis vectors in lab frame)
- `φj::Tφ`: Half-opening angle
- `z::Tz`: Axial extent interval (along-axis coordinate range)

# Constructor
Use `Conical(; axis, φj, z)` where `axis` is a unit vector defining the cone axis.
The rotation matrix is computed automatically.

# Matrix usage
- Lab → local: `r_local = R_local_to_lab' * r_lab` (transpose needed)
- Local → lab: `r_lab = R_local_to_lab * r_local` (direct)
"""
struct Conical{TR, Tφ, Tz} <: AbstractGeometry
	R_local_to_lab::TR
	φj::Tφ
	z::Tz

	# Inner constructor: validate that R_local_to_lab is a matrix, not a vector
	# This prevents old code from creating Conical with an axis instead of rotation matrix
	function Conical(R_local_to_lab::SMatrix{3,3}, φj, z)
		return new{typeof(R_local_to_lab), typeof(φj), typeof(z)}(R_local_to_lab, φj, z)
	end

	# Fallback for invalid types (like SVector) - convert axis to rotation matrix
	function Conical(axis::SVector{3}, φj, z)
		@assert isapprox(norm(axis), 1; atol=√eps(float(eltype(axis))))
		R = S._axis_to_rotation(axis)
		return new{typeof(R), typeof(φj), typeof(z)}(R, φj, z)
	end
end
Conical(; axis::AbstractVector, φj, z) = Conical(axis, φj, z)



"""
    InertialWorldline{TX, TU}

Inertial (constant-velocity) worldline: `x(τ) = x0 + u * τ`.

# Fields
- `x0::TX`: Reference spacetime event (`FourPosition`)
- `u::TU`: Constant 4-velocity (`FourVelocity`)
"""
@kwdef struct InertialWorldline{TX, TU}
	x0::TX
	u::TU
end

"""
    WorldtubeEllipsoid{TW, TS} <: AbstractGeometry

Ellipsoidal body moving along a worldline. The shape is defined in the rest frame.

# Fields
- `center::TW`: Worldline of the center of the ellipsoid (e.g., `InertialWorldline`)
- `sizes::TS`: Rest-frame semi-axes `SVector{3}` (x, y, z in rest frame)
"""
@kwdef struct Ellipsoid{TW, TS} <: AbstractGeometry
	center::TW
	sizes::TS
end

"""
    Cylindrical{TR, TRadius, Tz}

Cylindrical geometry with rotation matrix, fixed radius, and axial extent.

# Fields
- `R_local_to_lab::TR`: 3×3 rotation matrix (local jet frame → lab frame)
- `radius::TRadius`: Cylinder radius (constant)
- `z::Tz`: Axial extent interval

# Constructor
Use `Cylindrical(; axis, radius, z)` where `axis` is a unit vector.
"""
struct Cylindrical{TR, TRadius, Tz} <: AbstractGeometry
	R_local_to_lab::TR
	radius::TRadius
	z::Tz

	function Cylindrical(R_local_to_lab::SMatrix{3,3}, radius, z)
		return new{typeof(R_local_to_lab), typeof(radius), typeof(z)}(R_local_to_lab, radius, z)
	end

	function Cylindrical(axis::SVector{3}, radius, z)
		@assert isapprox(norm(axis), 1; atol=√eps(float(eltype(axis))))
		R = S._axis_to_rotation(axis)
		return new{typeof(R), typeof(radius), typeof(z)}(R, radius, z)
	end
end
Cylindrical(; axis::AbstractVector, radius, z) = Cylindrical(SVector{3}(axis), radius, z)

"""
    Ball{TW, TS} <: AbstractGeometry

Spherical body moving along a worldline. Same as `Ellipsoid` but with a single radius.
The `z_interval` implementation is significantly faster because it avoids the rest-frame
Lorentz boost, using Minkowski inner product invariance instead.

# Fields
- `center::TW`: Worldline of the center (e.g., `InertialWorldline`)
- `size::TS`: Rest-frame radius (scalar)
"""
@kwdef struct Ball{TW, TS} <: AbstractGeometry
	center::TW
	size::TS
end

end # module Geometries


"""Prepare for computations by caching trig values"""
prepare_for_computations(g::Geometries.Conical) = @modify(Geometries.AngleTrigCached_fromangle, g.φj)
ustrip(g::Geometries.Conical) = @modify(z -> _u_to_code(z, UCTX.L0), g.z)

prepare_for_computations(g::Geometries.Cylindrical) = g
ustrip(g::Geometries.Cylindrical) = @p let
	g
	@modify(_u_to_code(_, UCTX.L0), __.z)
	@modify(_u_to_code(_, UCTX.L0), __.radius)
end

prepare_for_computations(g::Geometries.Ellipsoid) = g

@inline four_velocity(g::Geometries.Ellipsoid, x4) = four_velocity(g.center)
@inline four_velocity(wl::Geometries.InertialWorldline) = wl.u

prepare_for_computations(g::Geometries.Ball) = g

@inline four_velocity(g::Geometries.Ball, x4) = four_velocity(g.center)

@inline function z_interval(g::Geometries.Ball, ray::Ray)
	# Worldtube of a rigid sphere moving with constant 4-velocity u.
	# For a sphere (all sizes equal), the quadratic can be solved directly from
	# Minkowski inner products, without Lorentz-boosting to the rest frame.
	# Key identity: for V ⊥ u (Minkowski), |v_rest|² = η(V,V).
	@assert g.center isa Geometries.InertialWorldline
	u = g.center.u
	Δ0 = ray.x0 - g.center.x0
	kdir = direction4(ray)

	a = minkowski_dot(u, kdir)
	b = minkowski_dot(u, Δ0)

	# η(P0,P0) = η(Δ0,Δ0) + b², η(P1,P1) = η(kdir,kdir) + a², η(P0,P1) = η(Δ0,kdir) + a·b
	p0_sq = minkowski_dot(Δ0, Δ0) + b^2
	p1_sq = minkowski_dot(kdir, kdir) + a^2
	p01   = minkowski_dot(Δ0, kdir) + a * b

	R_sq = g.size^2
	A = p1_sq / R_sq
	B = 2 * p01 / R_sq
	C = p0_sq / R_sq - one(A)
	D = B^2 - 4 * A * C

	if !(D > 0) || iszero(A)
		s0 = zero(A)
		return s0 .. (s0 - eps(oneunit(s0)))
	end

	sD = √(D)
	s1 = (-B - sD) / (2 * A)
	s2 = (-B + sD) / (2 * A)
	return min(s1, s2) .. max(s1, s2)
end

function z_interval(g::Geometries.Ellipsoid, ray::Ray)
	# Worldtube of a rigid axis-aligned ellipsoid moving with constant 4-velocity u.
	# In the comoving frame, define Δ⊥ as the displacement orthogonal to u, and require
	# (Δx/sx)^2 + (Δy/sy)^2 + (Δz/sz)^2 ≤ 1.
	# `sizes == SVector(R,R,R)` recovers the moving sphere.
	# Ray: x(s) = ray.x0 + s·d, where d = direction4(ray).
	# For SlowLight d = (1, n̂) (null); for FastLight d = (0, n̂) (spatial).
	@assert g.center isa Geometries.InertialWorldline
	u = g.center.u
	Δ0 = ray.x0 - g.center.x0
	kdir = direction4(ray)

	a = minkowski_dot(u, kdir)
	b = minkowski_dot(u, Δ0)
	P0 = Δ0 + u * b
	P1 = kdir + u * a

	p0 = _spatial_in_rest(u, P0)
	p1 = _spatial_in_rest(u, P1)

	s² = g.sizes .^ 2
	A = sum(p1 .^ 2 ./ s²)
	B = 2 * sum(p0 .* p1 ./ s²)
	C = sum(p0 .^ 2 ./ s²) - one(A)
	D = B^2 - 4 * A * C

	if !(D > 0) || iszero(A)
		s0 = zero(A)
		return s0 .. (s0 - eps(oneunit(s0)))
	end

	sD = √(D)
	s1 = (-B - sD) / (2 * A)
	s2 = (-B + sD) / (2 * A)
	return min(s1, s2) .. max(s1, s2)
end

"""Extract the axis vector (third column of rotation matrix) from a Conical geometry"""
@inline geometry_axis(g::Geometries.Conical) = g.R_local_to_lab[:,3]
Accessors.set(g::Geometries.Conical, ::typeof(geometry_axis), v::SVector{3}) = Geometries.Conical(v, g.φj, g.z)

@inline geometry_axis(g::Geometries.Cylindrical) = g.R_local_to_lab[:,3]
Accessors.set(g::Geometries.Cylindrical, ::typeof(geometry_axis), v::SVector{3}) = Geometries.Cylindrical(v, g.radius, g.z)

"""
Compute cylindrical coordinates from rotation matrix.

Takes rotation matrix R_local_to_lab and spacetime position x4.
Returns:
- `r`: spatial position in lab frame
- `z`: axial distance (z-coordinate in local frame)
- `r_local`: position in local frame (where z is along axis)
- `r_perp`: perpendicular component in lab frame (r - z*axis)
- `ρ`: cylindrical radius (distance from axis in transverse plane)
"""
@inline _cylindrical_coords(R, x4) = begin
	r = @swiz x4.xyz
	# if we needed r_local, that's probably optimal:
	# r_local = R' * r  # Transform to local frame (R' maps lab → local)
	# z = r_local[3]    # Axial coordinate is z-component in local frame
	# ρ = hypot(r_local[1], r_local[2])  # Cylindrical radius in transverse plane
	# # Perpendicular component: zero out z in local frame, transform back to lab
	# r_perp = R * SVector(r_local[1], r_local[2], zero(z))

	# but in practice, we don't need r_local:
	axis = R[:, 3]  # Third column is axis in lab frame
	z = dot(axis, r)
	r_perp = r - z * axis
	ρ = norm(r_perp)
	return (; r, z, r_perp, ρ)
end

"""
    natural_coords(g::Conical, x4) -> NamedTuple

Returns `(z, ρ, η)` where:
- `z`: axial coordinate (distance along axis)
- `ρ`: cylindrical radius (perpendicular distance from axis)
- `η`: normalized radial coordinate `η = ρ / (z * tan(φj))`
"""
@inline function natural_coords(g::Geometries.Conical, x4)
	(; z, ρ) = _cylindrical_coords(g.R_local_to_lab, x4)
	η = ρ / (z * tan(g.φj))
	return (; z, ρ, η)
end

"""
    natural_coords(g::Conical, x4, ::Val{:z}) -> Number

Returns just the axial coordinate `z` (optimized version).
"""
@inline natural_coords(g::Geometries.Conical, x4, ::Val{:z}) =
	dot(geometry_axis(g), @swiz x4.xyz)

"""
    is_inside(g::Conical, x4) -> Bool

Test if position is inside the conical volume.
"""
@inline function is_inside(g::Geometries.Conical, x4)
	(; z, ρ) = _cylindrical_coords(g.R_local_to_lab, x4)
	return (z ∈ g.z) && ρ ≤ z * tan(g.φj)
end

"""
    z_interval(g::Conical, ray) -> Interval

Compute ray-cone intersection interval in the ray parameter `s`.

The ray is parameterized as `r(s) = r0 + s·n̂` where `n̂ = direction3(ray)`.
Returns the interval of `s` values for which the ray is inside the truncated cone.
"""
function z_interval(g::Geometries.Conical, ray::Ray)
	@boundscheck @assert 0 ≤ leftendpoint(g.z) ≤ rightendpoint(g.z)
	c2 = cos(g.φj)^2
	r0 = @swiz ray.x0.xyz
	n̂ = direction3(ray)
	axis = geometry_axis(g)
	a_n = dot(axis, n̂)    # projection of ray direction onto cone axis

	# Cone inequality: (axis⋅r)² ≥ cos²(φj)·|r|².
	# With r(s) = r0 + s·n̂:
	#   axis⋅r(s) = α0 + a_n·s
	#   |r(s)|² = |r0|² + 2·d0·s + s²
	# where d0 = r0⋅n̂.
	α0 = dot(axis, r0)
	d0 = dot(r0, n̂)
	r02 = dot(r0, r0)
	A = a_n^2 - c2
	B = 2 * (α0 * a_n - c2 * d0)
	C = α0^2 - c2 * r02

	FT = eltype(ray.x0) |> float
	emptyseg = let
		sref = FT(0)
		sref .. (sref - eps(oneunit(sref)))
	end
	infseg = FT(-Inf) .. FT(Inf)

	# 1) Restrict to truncation interval in along-axis coordinate: axis⋅r(s) ∈ g.z.
	#    axis⋅r(s) = α0 + a_n·s, so s = (z_boundary - α0) / a_n.
	ss_trunc = if iszero(a_n)
		(α0 ∈ g.z) ? infseg : emptyseg
	else
		zmin, zmax = endpoints(g.z)
		s1 = (zmin - α0) / a_n
		s2 = (zmax - α0) / a_n
		(min(s1, s2) .. max(s1, s2))
	end

	isempty(ss_trunc) && return ss_trunc

	# 2) Intersect with the infinite cone.
	# Solve A·s² + B·s + C ≥ 0 and clip to `ss_trunc`.
	s_cone = let
		D = B^2 - 4 * A * C
		if iszero(A)
			# Linear case: B·s + C ≥ 0.
			if iszero(B)
				(C >= 0) ? infseg : emptyseg
			elseif B > 0
				(-C / B) .. infseg.right
			else
				infseg.left .. (-C / B)
			end
		elseif D < 0
			# No real roots: sign is constant (A and C have same sign when D<0).
			(C >= 0) ? infseg : emptyseg
		else
			s1 = (-B - √D) / (2 * A)
			s2 = (-B + √D) / (2 * A)
			slo, shi = minmax(s1, s2)
			if A < 0
				slo .. shi
			else
				# cone interior outside the roots
				left = infseg.left .. slo
				right = shi .. infseg.right
				# Prefer the side that overlaps `ss_trunc`
				L = ss_trunc ∩ left
				R = ss_trunc ∩ right
				@boundscheck @assert isempty(L) || isempty(R)
				!isempty(L) ? L : R
			end
		end
	end

	return ss_trunc ∩ s_cone
end


# ============================================================================
# Cylindrical geometry methods
# ============================================================================

@inline function natural_coords(g::Geometries.Cylindrical, x4)
	(; z, ρ) = _cylindrical_coords(g.R_local_to_lab, x4)
	η = ρ / g.radius
	return (; z, ρ, η)
end

@inline natural_coords(g::Geometries.Cylindrical, x4, ::Val{:z}) =
	dot(geometry_axis(g), @swiz x4.xyz)

@inline function is_inside(g::Geometries.Cylindrical, x4)
	(; z, ρ) = _cylindrical_coords(g.R_local_to_lab, x4)
	return (z ∈ g.z) && ρ ≤ g.radius
end

"""
    z_interval(g::Cylindrical, ray) -> Interval

Compute ray-cylinder intersection interval in the ray parameter `s`.

The ray is parameterized as `r(s) = r0 + s·n̂`.
Returns the interval of `s` values for which the ray is inside the truncated cylinder.
"""
function z_interval(g::Geometries.Cylindrical, ray::Ray)
	@boundscheck @assert 0 ≤ leftendpoint(g.z) ≤ rightendpoint(g.z)
	r0 = @swiz ray.x0.xyz
	n̂ = direction3(ray)
	axis = geometry_axis(g)
	a_n = dot(axis, n̂)

	FT = eltype(ray.x0) |> float
	emptyseg = let
		sref = FT(0)
		sref .. (sref - eps(oneunit(sref)))
	end
	infseg = FT(-Inf) .. FT(Inf)

	# 1) Axial truncation: axis⋅r(s) ∈ g.z
	α0 = dot(axis, r0)
	ss_trunc = if iszero(a_n)
		(α0 ∈ g.z) ? infseg : emptyseg
	else
		zmin, zmax = endpoints(g.z)
		s1 = (zmin - α0) / a_n
		s2 = (zmax - α0) / a_n
		(min(s1, s2) .. max(s1, s2))
	end

	isempty(ss_trunc) && return ss_trunc

	# 2) Cylinder: |r_perp(s)|² ≤ R²
	#    r_perp(s) = r0_perp + s·n̂_perp  (strip axial component)
	r0_perp = r0 - α0 * axis
	n̂_perp = n̂ - a_n * axis
	A = dot(n̂_perp, n̂_perp)
	B = 2 * dot(r0_perp, n̂_perp)
	C = dot(r0_perp, r0_perp) - g.radius^2

	ss_cyl = if iszero(A)
		# Ray parallel to axis: either always inside or always outside
		(C ≤ 0) ? infseg : emptyseg
	else
		D = B^2 - 4 * A * C
		if !(D > 0)
			emptyseg
		else
			sD = √(D)
			s1 = (-B - sD) / (2 * A)
			s2 = (-B + sD) / (2 * A)
			min(s1, s2) .. max(s1, s2)
		end
	end

	return ss_trunc ∩ ss_cyl
end

rotation_local_to_lab(geom::Geometries.Cylindrical) = geom.R_local_to_lab

# ============================================================================
# Coordinate transformation helpers
# ============================================================================

"""Helper: minimal rotation from unit vector a to unit vector b"""
@inline _rotate_minimal(a::SVector{3}, b::SVector{3}, v::SVector{3}) = begin
	ET = float(promote_type(eltype(a), eltype(b), eltype(v)))

	c = dot(a, b)
	vx = cross(a, b)
	s = norm(vx)

	# Rodrigues rotation formula for the minimal rotation mapping â -> b̂.
	# Special-case (anti)parallel vectors for numerical stability.
	if s < √eps(s)
		if c > 0
			return ET.(v)
		else
			# 180° rotation: axis is not unique; pick a deterministic one perpendicular to â.
			axis = abs(a.y) < 0.9 ? SVector{3,ET}(0, 1, 0) : SVector{3,ET}(1, 0, 0)
			k = normalize(cross(a, axis))
			# For θ = π: R(v) = v - 2 k×(k×v) = 2(k⋅v)k - v.
			return 2 * dot(k, v) * k - v
		end
	end

	k = vx / s
	# v_rot = v cosθ + (k×v) sinθ + k(k⋅v)(1-cosθ)
	return v * c + cross(k, v) * s + k * dot(k, v) * (1 - c)
end

"""
    rotation_local_to_lab(geom) -> SMatrix{3,3}

Return the 3×3 rotation matrix `R` that maps local jet frame to lab frame.

Matrix columns are lab-frame images of the local basis vectors.
Constructed via minimal rotation from lab `ẑ` axis to `geometry_axis(geom)`.

Usage:
- lab → local coordinates: `r_local = R' * r_lab` (transpose needed)
- local → lab coordinates: `r_lab = R * r_local` (direct)
"""
rotation_local_to_lab(geom::Geometries.Conical) = geom.R_local_to_lab

function _axis_to_rotation(axis::SVector{3})
	ẑ = SVector(0, 0, 1)
	x̂ = SVector(1, 0, 0)
	ŷ = SVector(0, 1, 0)

	ex = _rotate_minimal(ẑ, axis, x̂)
	ey = _rotate_minimal(ẑ, axis, ŷ)
	ez = axis

	return hcat(ex, ey, ez)
end

"""
    rotate_lab_to_local(geom, r::SVector{3}) -> SVector{3}

Convert a lab-frame spatial 3-position `r` to local coordinates defined by
the geometry rotation matrix.
"""
@inline rotate_lab_to_local(geom::Geometries.AbstractGeometry, r::SVector{3}) = rotation_local_to_lab(geom)' * r

"""
    rotate_local_to_lab(geom, r_local::SVector{3}) -> SVector{3}

Convert a local-frame spatial 3-position back to lab coordinates.
"""
@inline rotate_local_to_lab(geom::Geometries.AbstractGeometry, r_local::SVector{3}) = rotation_local_to_lab(geom) * r_local

"""
	ray_in_local_coords(ray, geom; s_range) -> StructArray{SVector{3}}

Compute the two endpoints of a ray segment in the geometry's local coordinate frame.

# Arguments
- `ray`: a `Ray`.
- `geom`: geometry (or `EmissionRegion`) whose local frame is used.
- `s_range`: `Interval` of the ray parameter `s` (distance along the ray direction).

# Returns
A 2-element `StructArray` of `SVector{3}`: the local-frame positions at the two `s_range` endpoints.
The local z-axis is aligned with `geometry_axis(geom)`.
"""
function ray_in_local_coords(ray::Ray, geom; s_range)
	n̂ = direction3(ray)
	r0 = @swiz ray.x0.xyz
	map(endpoints(s_range)) do s
		r_lab = r0 + s * n̂
		rotate_lab_to_local(geom, r_lab)
	end |> collect |> StructArray
end

"""
	ray_in_local_coords(ray::RayGR2, geom; s_range) -> StructArray{SVector{3}}

Compute the GR-deflected ray path in the geometry's local coordinate frame.

Returns 4 points tracing the bent photon path (or 2 if the photon is captured):
1. Far observer end of the incoming ray
2. Near-BH end of the incoming ray
3. Near-BH end of the outgoing ray
4. Far source end of the outgoing ray
"""
function ray_in_local_coords(ray::RayGR2, geom; s_range)
	(; ray_in, ray_out, bh_position) = ray
	s_lo, s_hi = extrema(endpoints(s_range))

	n̂_in = direction3(ray_in)
	r0_in = @swiz ray_in.x0.xyz
	s_bh_in = _s_closest_to_point(ray_in, bh_position)

	pt_obs   = r0_in + s_hi * n̂_in      # far observer end (n̂_in points toward observer)
	pt_bh_in = r0_in + s_bh_in * n̂_in   # near BH on incoming ray

	if is_captured_ray(ray_out)
		return map(SVector(pt_obs, pt_bh_in)) do p
			rotate_lab_to_local(geom, p)
		end |> collect |> StructArray
	end

	n̂_out = direction3(ray_out)
	r0_out = @swiz ray_out.x0.xyz
	s_bh_out = _s_closest_to_point(ray_out, bh_position)

	pt_bh_out = r0_out + s_bh_out * n̂_out  # near BH on outgoing ray
	pt_src    = r0_out + s_lo * n̂_out       # far source end (n̂_out points toward BH, so negative s → source)

	map(SVector(pt_obs, pt_bh_in, pt_bh_out, pt_src)) do p
		rotate_lab_to_local(geom, p)
	end |> collect |> StructArray
end

"""
	camera_fov_in_local_coords(cam, geom; v=0, s_range) -> StructArray{SVector{3}}

Compute the four corners of the camera's field-of-view rectangle in the geometry's local frame.

The rectangle is the cross-section of the camera ray bundle at a fixed screen-plane offset `v`
(along `cam.e2`), swept over the full `u`-range of `cam.xys` and the full `s_range` depth.

# Arguments
- `cam`: a `Camera`.
- `geom`: geometry (or `EmissionRegion`) whose local frame is used.
- `v`: fixed offset along `cam.e2` (default `0`).
- `s_range`: `Interval` of the ray parameter `s` (distance along `cam.n`).

# Returns
A 4-element `StructArray` of `SVector{3}`: the local-frame positions of the four corners,
ordered `(umin,s1), (umax,s1), (umax,s2), (umin,s2)`.
"""
@unstable function camera_fov_in_local_coords(cam, geom; v=0, s_range)
	umin, umax = extrema(uv -> uv[1], cam.xys)
	s1, s2 = endpoints(s_range)

	corners_lab = map([(umin, s1), (umax, s1), (umax, s2), (umin, s2)]) do (u, s)
		cam.origin + u * cam.e1 + v * cam.e2 + s * cam.n
	end |> StructArray

	map(corners_lab) do r_lab
		rotate_lab_to_local(geom, r_lab)
	end
end

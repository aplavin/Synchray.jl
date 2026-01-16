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

end # module Geometries


"""Prepare for computations by caching trig values"""
prepare_for_computations(g::Geometries.Conical) = @modify(Geometries.AngleTrigCached_fromangle, g.φj)
ustrip(g::Geometries.Conical) = @modify(z -> _u_to_code(z, UCTX.L0), g.z)

"""Extract the axis vector (third column of rotation matrix) from a Conical geometry"""
@inline geometry_axis(g::Geometries.Conical) = g.R_local_to_lab[:,3]
Accessors.set(g::Geometries.Conical, ::typeof(geometry_axis), v::SVector{3}) = Geometries.Conical(v, g.φj, g.z)

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

Compute ray-cone intersection interval in lab `z` coordinate.
"""
function z_interval(g::Geometries.Conical, ray)
	# RayZ is a line of sight parameterized by z: r(z) = r0 + z*e_z.
	# Assumes standard camera convention: ray.x0.z == 0 and cone apex at origin.
	@assert iszero(ray.x0.z)
	@assert 0 ≤ leftendpoint(g.z) ≤ rightendpoint(g.z)
	
	c2 = cos(g.φj)^2
	r0 = SVector(ray.x0.x, ray.x0.y, ray.x0.z)
	axis = geometry_axis(g)
	az = axis.z

	# Cone inequality: (a⋅r)^2 ≥ c^2 * (r⋅r).
	# With r(z)=r0+z*e_z and r0⋅e_z=0:
	# (α0 + az*z)^2 - c^2*(|r0|^2 + z^2) >= 0
	α0 = dot(axis, r0)
	r02 = dot(r0, r0)
	A = az^2 - c2
	B = 2 * α0 * az
	C = α0^2 - c2 * r02

	FT = eltype(ray.x0) |> float
	emptyseg = let
		zref = FT(ray.x0.z)
		zref .. (zref - eps(zref))
	end
	infseg = FT(-Inf) .. FT(Inf)

	# 1) Restrict to truncation interval in along-axis coordinate: z_axis(z_lab) ∈ g.z.
	zs = if iszero(az)
		(α0 ∈ g.z) ? infseg : emptyseg
	else
		zmin, zmax = endpoints(g.z)
		z1 = (zmin - α0) / az
		z2 = (zmax - α0) / az
		(min(z1, z2) .. max(z1, z2))
	end

	isempty(zs) && return zs

	# 2) Intersect with the infinite cone (already symmetric); half-cone handled above.
	# Solve A z^2 + B z + C ≥ 0 and clip to `zs`.
	z_cone = let
		D = B^2 - 4 * A * C
		if iszero(A)
			# Linear case: B z + C ≥ 0.
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
			z1 = (-B - √D) / (2 * A)
			z2 = (-B + √D) / (2 * A)
			zlo, zhi = minmax(z1, z2)
			if A < 0
				zlo .. zhi
			else
				# cone interior outside the roots
				@assert zlo ≤ 0 ≤ zhi  (zlo, zhi)
				left = infseg.left .. zlo
				right = zhi .. infseg.right
				# Prefer the side that overlaps `zs`
				L = zs ∩ left
				R = zs ∩ right
				@assert isempty(L) || isempty(R)
				!isempty(L) ? L : R
			end
		end
	end

	return zs ∩ z_cone
end

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
    ray_in_local_coords(ray, geom; z_range) -> StructArray{SVector{3}}

Project a ray into local coordinates.

Returns 3-vectors in local frame defining the line segment for the given `z_range` interval.
Use `@swiz` at call sites to extract desired components, e.g., `@swiz r.zx` for (s, h) coordinates.
"""
function ray_in_local_coords(ray, geom; z_range)
	map(endpoints(z_range)) do z
		r_lab = SVector(ray.x0.x, ray.x0.y, z)
		rotate_lab_to_local(geom, r_lab)
	end |> collect |> StructArray
end

"""
    camera_fov_in_local_coords(cam, geom; y=0, z_range) -> StructArray{SVector{3}}

Project the camera's field of view into local coordinates.

Returns four corners of the band (polygon) swept by camera rays in the y=y plane.
Use `@swiz` at call sites to extract desired components.
"""
@unstable function camera_fov_in_local_coords(cam, geom; y=0, z_range)
	xmin, xmax = extrema(xy -> xy.x, cam.xys)
	z1, z2 = endpoints(z_range)

	corners_lab = [
		SVector(xmin, y, z1),
		SVector(xmax, y, z1),
		SVector(xmax, y, z2),
		SVector(xmin, y, z2),
	] |> StructArray

	map(corners_lab) do r_lab
		rotate_lab_to_local(geom, r_lab)
	end
end

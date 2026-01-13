"""
Module for geometry types and coordinate transformations.

Geometries define the shape, boundaries, and natural coordinate systems
for emission regions.
"""
module Geometries

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
	tan::T
	cos::T
end

Base.tan(x::AngleTrigCached) = x.tan
Base.cos(x::AngleTrigCached) = x.cos
AngleTrigCached_fromangle(φ) = AngleTrigCached(tan(φ), cos(φ))

"""
    Conical{Ta, Tφ, Tz}

Conical geometry with axis, half-opening angle, and axial extent.

# Fields
- `axis::Ta`: Unit vector defining the cone axis
- `φj::Tφ`: Half-opening angle
- `z::Tz`: Axial extent interval (along-axis coordinate range)
"""
@kwdef struct Conical{Ta, Tφ, Tz} <: AbstractGeometry
	axis::Ta
	φj::Tφ
	z::Tz
end

end # module Geometries


"""Prepare for computations by caching trig values"""
prepare_for_computations(g::Geometries.Conical) = Geometries.Conical(g.axis, Geometries.AngleTrigCached_fromangle(g.φj), g.z)

@accessor geometry_axis(g::Geometries.Conical) = g.axis

"""
Compute cylindrical coordinates relative to axis.
Returns: spatial position `r`, axial distance `z`, perpendicular component `r_perp`, cylindrical radius `ρ`.
"""
@inline _cylindrical_coords(axis, x4) = begin
	r = SVector(x4.x, x4.y, x4.z)
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
function natural_coords(g::Geometries.Conical, x4)
	(; z, ρ) = _cylindrical_coords(g.axis, x4)
	η = ρ / (z * tan(g.φj))
	return (; z, ρ, η)
end

"""
    natural_coords(g::Conical, x4, ::Val{:z}) -> Number

Returns just the axial coordinate `z` (optimized version).
"""
function natural_coords(g::Geometries.Conical, x4, ::Val{:z})
	r = SVector(x4.x, x4.y, x4.z)
	return dot(g.axis, r)
end

"""
    is_inside(g::Conical, x4) -> Bool

Test if position is inside the conical volume.
"""
function is_inside(g::Geometries.Conical, x4)
	(; z, ρ) = _cylindrical_coords(g.axis, x4)
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
	az = g.axis.z

	# Cone inequality: (a⋅r)^2 ≥ c^2 * (r⋅r).
	# With r(z)=r0+z*e_z and r0⋅e_z=0:
	# (α0 + az*z)^2 - c^2*(|r0|^2 + z^2) >= 0
	α0 = dot(g.axis, r0)
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
    rotation_lab_to_local(geom) -> SMatrix{3,3}

Return the 3×3 rotation matrix `R` for minimal rotation from lab frame that maps
lab `ẑ` axis to `geometry_axis(geom)`.

Matrix columns are lab-frame images of the local basis vectors.

- lab → local coordinates: `r_local = R' * r_lab`
- local → lab coordinates: `r_lab = R * r_local`
"""
function rotation_lab_to_local(geom::Geometries.AbstractGeometry)
	axis = geometry_axis(geom)
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
rotate_lab_to_local(geom::Geometries.AbstractGeometry, r::SVector{3}) = 
	rotation_lab_to_local(geom)' * r

"""
    rotate_local_to_lab(geom, r_local::SVector{3}) -> SVector{3}

Convert a local-frame spatial 3-position back to lab coordinates.
"""
rotate_local_to_lab(geom::Geometries.AbstractGeometry, r_local::SVector{3}) = 
	rotation_lab_to_local(geom) * r_local

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

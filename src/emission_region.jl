"""
    EmissionRegion{G, Tu, Temission} <: AbstractMedium

Generic emission region combining geometry, velocity field, and emission model.

# Fields
- `geometry::G`: Geometry defining shape and coordinates (e.g. `Conical`, `WorldtubeEllipsoid`)
- `velocity::Tu`: Velocity specification (`VelocitySpec`, `CombinedVelocity`, or `nothing` for worldline-based geometries)
- `emission::Temission`: Emission model (e.g. `SynchrotronEmission`)

# Examples
```julia
# Jet with explicit velocity field:
EmissionRegion(
    geometry = Geometries.Conical(axis = SVector(0, 0, 1), φj = 0.1, z = 1.0 .. 10.0),
    velocity = VelocitySpec(Directions.Axial(), Profiles.Constant(10.0)),
    emission = SynchrotronEmission(
        ne = Profiles.Axial(PowerLaw(-2; val0=1e4, s0=1.0)),
        B = BFieldSpec(Profiles.Axial(PowerLaw(-1; val0=1e-3, s0=1.0)), Directions.Scalar(), b -> FullyTangled(b)),
        electrons = IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=0.1)
    )
)

# Moving blob (velocity derived from worldline):
EmissionRegion(
    geometry = Geometries.WorldtubeEllipsoid(
        worldline = Geometries.InertialWorldline(FourPosition(0,0,0,5), construct(FourVelocity, gamma => 10, direction3 => SA[0,0,1])),
        sizes = SA[0.5, 0.5, 0.5]
    ),
    emission = SynchrotronEmission(
        ne = Profiles.Constant(1e4),
        B = BFieldSpec(Profiles.Constant(1e-3), Directions.Scalar(), b -> FullyTangled(b)),
        electrons = IsotropicPowerLawElectrons(; p=2.5)
    )
)
```
"""
@kwdef struct EmissionRegion{G, Tu, Temission} <: AbstractMedium
    geometry::G
    velocity::Tu = nothing
    emission::Temission
end

# Delegate geometry methods
@inline z_interval(region::EmissionRegion, ray) = z_interval(region.geometry, ray)
@inline is_inside(region::EmissionRegion, x4) = is_inside(region.geometry, x4)
@unstable @accessor geometry_axis(region::EmissionRegion) = geometry_axis(region.geometry)

# Velocity dispatch: explicit field or derived from worldline
@inline four_velocity(region::EmissionRegion, x4) = four_velocity(region.velocity, region.geometry, x4)
@inline four_velocity(::Nothing, geom, x4) = four_velocity(geom, x4)

# For VelocitySpec: direction + magnitude from profile
@inline function four_velocity(vel::VelocitySpec, geom, x4)
    speed_val = vel.scale(geom, x4)
    v_dir = field_direction(vel.dir, geom, x4)
    return construct(FourVelocity, vel.kind => speed_val, direction3 => v_dir)
end

# For CombinedVelocity: add β vectors in lab frame
@inline function four_velocity(vel::CombinedVelocity, geom, x4)
    u1 = four_velocity(vel.v1, geom, x4)
    u2 = four_velocity(vel.v2, geom, x4)
    βv1 = beta(u1)
    βv2 = beta(u2)
    βv_total = βv1 + βv2
    return FourVelocity(βv_total)
end

# Delegate to emission model
@inline emissivity_absorption(region::EmissionRegion, x4, k′) = emissivity_absorption(region.emission, region.geometry, x4, k′)
@inline emissivity_absorption_polarized(region::EmissionRegion, x4, k′) = emissivity_absorption_polarized(region.emission, region.geometry, x4, k′)

# Convenience accessors (delegate through emission)
@inline electron_density(region::EmissionRegion, x4) = electron_density(region.emission, region.geometry, x4)
@inline magnetic_field(region::EmissionRegion, x4) = magnetic_field(region.emission, region.geometry, x4)

# prepare_for_computations and ustrip propagation
@unstable prepare_for_computations(region::EmissionRegion) =
    modify(prepare_for_computations, region, @o _.geometry _.velocity _.emission)

@unstable ustrip(region::EmissionRegion) =
    modify(ustrip, region, @o _.geometry _.velocity _.emission)

# Visualization helpers forward to geometry
rotation_local_to_lab(region::EmissionRegion) = rotation_local_to_lab(region.geometry)
rotate_lab_to_local(region::EmissionRegion, r::SVector{3}) = rotate_lab_to_local(region.geometry, r)
rotate_local_to_lab(region::EmissionRegion, r_local::SVector{3}) = rotate_local_to_lab(region.geometry, r_local)
ray_in_local_coords(ray, region::EmissionRegion; s_range) = ray_in_local_coords(ray, region.geometry; s_range)
@unstable camera_fov_in_local_coords(cam, region::EmissionRegion; v=0, s_range) = camera_fov_in_local_coords(cam, region.geometry; v, s_range)

"""
    EmissionRegion{G, Tne, TB, Tu, Tmodel} <: AbstractSynchrotronMedium

Generic emission region with pluggable geometry and field profiles.

# Fields
- `geometry::G`: Geometry defining shape and coordinates
- `ne::Tne`: Electron density profile, callable as `ne(geom, x4) -> scalar`
- `B::TB`: Magnetic field specification (`BFieldSpec`)
- `velocity::Tu`: Velocity field specification (`VelocitySpec`)
- `model::Tmodel`: Synchrotron electron model

# Example
```julia
region = EmissionRegion(
    geometry = Geometries.Conical(
        axis = SVector(0, 0, 1),
        φj = 0.1,
        z = 1.0 .. 10.0
    ),
    ne = Profiles.Axial(PowerLaw(-2; val0=1e4, s0=1.0)),
    B = BFieldSpec(
        Profiles.Axial(PowerLaw(-1; val0=1e-3, s0=1.0)),
        Directions.Scalar(),
        b -> FullyTangled(b)
    ),
    velocity = VelocitySpec(
        Directions.Axial(),
        Profiles.Constant(10.0)
    ),
    model = IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=0.1)
)
```
"""
@kwdef struct EmissionRegion{G, Tne, TB, Tu, Tmodel} <: AbstractSynchrotronMedium
    geometry::G
    ne::Tne
    B::TB
    velocity::Tu
    model::Tmodel
end

# Delegate geometry methods
@inline z_interval(region::EmissionRegion, ray) = z_interval(region.geometry, ray)
@inline is_inside(region::EmissionRegion, x4) = is_inside(region.geometry, x4)
@unstable @accessor geometry_axis(region::EmissionRegion) = geometry_axis(region.geometry)

# Implement medium interface
@inline electron_density(region::EmissionRegion, x4) = region.ne(region.geometry, x4)

@inline function magnetic_field(region::EmissionRegion, x4)
    g = region.geometry
    b_mag = region.B.scale(g, x4)
    b_dir = field_direction(region.B.dir, g, x4)
    return region.B.wrap(b_mag * b_dir)
end

@inline function four_velocity(region::EmissionRegion, x4)
    g = region.geometry
    speed_val = region.velocity.scale(g, x4)
    v_dir = field_direction(region.velocity.dir, g, x4)
    return construct(FourVelocity, region.velocity.kind => speed_val, direction => v_dir)
end

@inline synchrotron_model(region::EmissionRegion) = region.model

# prepare_for_computations propagation
@unstable prepare_for_computations(region::EmissionRegion) = modify(prepare_for_computations, region, @o _.geometry _.ne _.B _.velocity _.model)
@unstable ustrip(region::EmissionRegion) = @p let
    region
    modify(ustrip, __, @o _.geometry _.B _.velocity _.model)
    @modify(ustrip(_; valu=UCTX.ne0), __.ne)
end

# Visualization helpers forward to geometry
rotation_lab_to_local(region::EmissionRegion) = rotation_lab_to_local(region.geometry)
rotate_lab_to_local(region::EmissionRegion, r::SVector{3}) = rotate_lab_to_local(region.geometry, r)
rotate_local_to_lab(region::EmissionRegion, r_local::SVector{3}) = rotate_local_to_lab(region.geometry, r_local)
ray_in_local_coords(ray, region::EmissionRegion; z_range) = ray_in_local_coords(ray, region.geometry; z_range)
@unstable camera_fov_in_local_coords(cam, region::EmissionRegion; y=0, z_range) = camera_fov_in_local_coords(cam, region.geometry; y, z_range)

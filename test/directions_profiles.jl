@testitem "Directions module" begin
    import Synchray as S

    @test Directions.Scalar() isa Directions.AbstractDirection
    @test Directions.Axial() isa Directions.AbstractDirection
    @test Directions.Radial() isa Directions.AbstractDirection
    @test Directions.Toroidal() isa Directions.AbstractDirection
    @test Directions.Helical(π/4) isa Directions.AbstractDirection
    
    # Scalar returns 1
    geom = nothing  # not used for Scalar
    x4 = S.FourPosition(0, 0, 0, 1)
    @test S.field_direction(Directions.Scalar(), geom, x4) == 1
end

@testitem "Profiles module" begin
    import Synchray as S
    # Test profile types exist and are callable
    @test Profiles.Axial(s -> s^-2) isa Profiles.Axial
    @test Profiles.Transverse(η -> exp(-η^2)) isa Profiles.Transverse
    @test Profiles.AxialTransverse(s -> s^-2, η -> exp(-η^2)) isa Profiles.AxialTransverse
    @test Profiles.Natural(c -> c.z * c.η) isa Profiles.Natural
    @test Profiles.Raw((g, x4) -> 1.0) isa Profiles.Raw
    @test Profiles.Constant(42.0) isa Profiles.Constant
    @test Profiles.Modified(Profiles.Constant(1.0), (g, x4, v) -> 2v) isa Profiles.Modified
end

@testitem "Profile evaluation with geometry" begin
    import Synchray as S
    using Synchray.Geometries

    geom = Geometries.Conical(SVector(0, 0, 1), 0.1, 1.0 .. 5.0)
    x4 = S.FourPosition(0, 0.05, 0, 2.0)
    
    # Constant profile
    p_const = Profiles.Constant(42.0)
    @test p_const(geom, x4) == 42.0
    
    # Axial profile
    p_axial = Profiles.Axial(s -> s^2)
    @test p_axial(geom, x4) ≈ 4.0
    
    # Transverse profile
    coords = S.natural_coords(geom, x4)
    p_trans = Profiles.Transverse(η -> 2η)
    @test p_trans(geom, x4) ≈ 2 * coords.η
    
    # AxialTransverse profile
    p_prod = Profiles.AxialTransverse(s -> s^2, η -> 3η)
    @test p_prod(geom, x4) ≈ 4.0 * (3 * coords.η)
    
    # Natural profile
    p_nat = Profiles.Natural(c -> c.z + c.η)
    @test p_nat(geom, x4) ≈ coords.z + coords.η
    
    # Raw profile
    p_raw = Profiles.Raw((g, x4) -> dot(S.geometry_axis(g), SVector(x4.x, x4.y, x4.z)))
    @test p_raw(geom, x4) ≈ 2.0
    
    # Modified profile
    p_base = Profiles.Constant(10.0)
    p_mod = Profiles.Modified(p_base, (g, x4, v) -> 2v)
    @test p_mod(geom, x4) == 20.0
end

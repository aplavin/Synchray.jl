@testitem "Directions module" begin
    import Synchray as S

    @test Directions.Scalar() isa Directions.AbstractDirection
    @test Directions.Axial() isa Directions.AbstractDirection
    @test Directions.Radial() isa Directions.AbstractDirection
    @test Directions.Toroidal() isa Directions.AbstractDirection
    @test Directions.HelicalAT(π/4) isa Directions.AbstractDirection
    @test Directions.HelicalRT(π/4) isa Directions.AbstractDirection

    # Scalar returns 1
    geom = nothing  # not used for Scalar
    x4 = S.FourPosition(0, 0, 0, 1)
    @test S.field_direction(Directions.Scalar(), geom, x4) == 1
end

@testitem "HelicalRT direction" begin
    import Synchray as S

    # Setup: conical geometry with z-axis
    geom = S.Geometries.Conical(SVector(0, 0, 1), 0.1, 1.0 .. 10.0)
    axis = SVector(0.0, 0.0, 1.0)

    # Test position off-axis
    x4 = S.FourPosition(0.0, 1.0, 0.0, 5.0)  # at (x=1, y=0, z=5)

    # Expected basis vectors at this position
    e_r_expected = SVector(1.0, 0.0, 0.0)  # radial: pointing in +x direction
    e_phi_expected = SVector(0.0, 1.0, 0.0)  # toroidal: pointing in +y direction

    @testset "ψ=0 gives pure radial" begin
        dir = Directions.HelicalRT(0.0)
        v = S.field_direction(dir, geom, x4)
        @test norm(v) ≈ 1.0
        @test v ≈ e_r_expected atol=1e-10
    end

    @testset "ψ=π/2 gives pure toroidal" begin
        dir = Directions.HelicalRT(π/2)
        v = S.field_direction(dir, geom, x4)
        @test norm(v) ≈ 1.0
        @test v ≈ e_phi_expected atol=1e-10
    end

    @testset "ψ=π/4 gives equal mix" begin
        dir = Directions.HelicalRT(π/4)
        v = S.field_direction(dir, geom, x4)
        @test norm(v) ≈ 1.0
        expected = (e_r_expected + e_phi_expected) / sqrt(2)
        @test v ≈ expected atol=1e-10
    end

    @testset "Result is perpendicular to axis" begin
        for ψ in [0.0, π/6, π/4, π/3, π/2]
            dir = Directions.HelicalRT(ψ)
            v = S.field_direction(dir, geom, x4)
            @test abs(dot(v, axis)) < 1e-10
        end
    end

    @testset "On-axis returns zero" begin
        x4_on_axis = S.FourPosition(0.0, 0.0, 0.0, 5.0)
        dir = Directions.HelicalRT(π/4)
        v = S.field_direction(dir, geom, x4_on_axis)
        @test norm(v) ≈ 0.0
    end
end

@testitem "Profiles module" begin
    import Synchray as S
    # Test profile types exist and are callable
    @test Profiles.Axial(s -> s^-2) isa Profiles.Axial
    @test Profiles.Transverse(η -> exp(-η^2)) isa Profiles.Transverse
    @test Profiles.Radial(r -> exp(-r^2)) isa Profiles.Radial
    @test Profiles.AxialTransverse(s -> s^-2, η -> exp(-η^2)) isa Profiles.AxialTransverse
    @test Profiles.Natural(c -> c.z * c.η) isa Profiles.Natural
    @test Profiles.Raw((g, x4) -> 1.0) isa Profiles.Raw
    @test Profiles.Constant(42.0) isa Profiles.Constant
    @test Profiles.Modified(Profiles.Constant(1.0), (g, x4, v) -> 2v) isa Profiles.Modified
    @test Profiles.LinearInterp(((1.0, 10.0), (2.0, 20.0))) isa Profiles.LinearInterp
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
    
    # Radial profile
    r = norm(@swiz x4.xyz)
    p_radial = Profiles.Radial(r -> 3r)
    @test p_radial(geom, x4) ≈ 3 * r
    
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

@testitem "LinearInterp profile" begin
    import Synchray as S
    
    # Basic construction
    li = Profiles.LinearInterp(((1.0, 10.0), (3.0, 30.0), (5.0, 50.0)))
    @test li isa Profiles.LinearInterp
    
    # Auto-sorting: provide unsorted points
    li_unsorted = Profiles.LinearInterp(((5.0, 50.0), (1.0, 10.0), (3.0, 30.0)))
    @test li_unsorted == li  # Should be equal after sorting
    
    # Flat value before first point
    @test li(0.5) == 10.0
    @test li(-1.0) == 10.0
    
    # Exact values at points
    @test li(1.0) == 10.0
    @test li(3.0) == 30.0
    @test li(5.0) == 50.0
    
    # Linear interpolation between points
    @test li(2.0) == 20.0  # Midpoint between 1.0 and 3.0
    @test li(1.5) == 15.0  # Quarter way
    @test li(4.0) == 40.0  # Midpoint between 3.0 and 5.0
    
    # Flat value after last point
    @test li(6.0) == 50.0
    @test li(10.0) == 50.0
    
    # Single point case
    li_single = Profiles.LinearInterp(((2.0, 100.0),))
    @test li_single(0.0) == 100.0
    @test li_single(2.0) == 100.0
    @test li_single(5.0) == 100.0
    
    # Two points case
    li_two = Profiles.LinearInterp(((0.0, 0.0), (10.0, 100.0)))
    @test li_two(-1.0) == 0.0
    @test li_two(0.0) == 0.0
    @test li_two(5.0) == 50.0
    @test li_two(10.0) == 100.0
    @test li_two(15.0) == 100.0
end

@testitem "LinearInterp with Axial profile" begin
    import Synchray as S
    using Synchray.Geometries
    
    # Create geometry and position
    geom = Geometries.Conical(SVector(0, 0, 1), 0.1, 1.0 .. 5.0)
    
    # Create LinearInterp and use with Axial
    li = Profiles.LinearInterp(((1.0, 100.0), (3.0, 200.0), (5.0, 150.0)))
    p_axial = Profiles.Axial(li)
    
    # Test evaluation at different positions
    x4_at_2 = S.FourPosition(0, 0.05, 0, 2.0)
    @test p_axial(geom, x4_at_2) ≈ 150.0  # Interpolated between 100 and 200
    
    x4_at_4 = S.FourPosition(0, 0.05, 0, 4.0)
    @test p_axial(geom, x4_at_4) ≈ 175.0  # Interpolated between 200 and 150
end

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

    # Test position off-axis: (x=1, y=0, z=5)
    x4 = S.FourPosition(0.0, 1.0, 0.0, 5.0)
    r = SVector(1.0, 0.0, 5.0)

    # Expected basis vectors at this position
    # Radial: from origin (NOT from axis!)
    e_r_expected = r / norm(r)  # ≈ (0.196, 0, 0.981)
    # Toroidal: azimuthal around axis
    e_phi_expected = SVector(0.0, 1.0, 0.0)

    @testset "ψ=0 gives pure radial (from origin)" begin
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

    @testset "ψ=π/4 gives normalized mix" begin
        dir = Directions.HelicalRT(π/4)
        v = S.field_direction(dir, geom, x4)
        @test norm(v) ≈ 1.0
        # e_r and e_phi are NOT orthogonal, so must normalize
        raw = (e_r_expected + e_phi_expected) / sqrt(2)
        expected = raw / norm(raw)
        @test v ≈ expected atol=1e-10
    end

    @testset "Result is always unit vector" begin
        for ψ in [0.0, π/6, π/4, π/3, π/2]
            dir = Directions.HelicalRT(ψ)
            v = S.field_direction(dir, geom, x4)
            @test norm(v) ≈ 1.0
        end
    end

    @testset "At origin returns zero" begin
        x4_origin = S.FourPosition(0.0, 0.0, 0.0, 0.0)
        dir = Directions.HelicalRT(π/4)
        v = S.field_direction(dir, geom, x4_origin)
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

@testitem "RigidRotation profile" begin
    import Synchray as S
    using Synchray.Geometries

    # Conical geometry with z-axis
    geom = Geometries.Conical(SVector(0, 0, 1), 0.1, 1.0 .. 10.0)

    # RigidRotation: β = β_ref × (ρ / ρ_ref)
    β_ref = 0.1
    ρ_ref = 2.0
    p = Profiles.RigidRotation(β_ref, ρ_ref)

    @test p isa Profiles.RigidRotation

    @testset "β scales linearly with ρ" begin
        # At ρ = ρ_ref: β = β_ref
        x4_at_ref = S.FourPosition(0.0, ρ_ref, 0.0, 5.0)  # ρ = 2.0
        @test p(geom, x4_at_ref) ≈ β_ref

        # At ρ = 2 × ρ_ref: β = 2 × β_ref
        x4_at_2ref = S.FourPosition(0.0, 2*ρ_ref, 0.0, 5.0)  # ρ = 4.0
        @test p(geom, x4_at_2ref) ≈ 2 * β_ref

        # At ρ = 0.5 × ρ_ref: β = 0.5 × β_ref
        x4_at_half = S.FourPosition(0.0, 0.5*ρ_ref, 0.0, 5.0)  # ρ = 1.0
        @test p(geom, x4_at_half) ≈ 0.5 * β_ref
    end

    @testset "At ρ = 0: β = 0" begin
        x4_on_axis = S.FourPosition(0.0, 0.0, 0.0, 5.0)
        @test p(geom, x4_on_axis) ≈ 0.0
    end

    @testset "Different z, same ρ gives same result" begin
        # ρ depends only on distance from axis, not on z
        x4_z5 = S.FourPosition(0.0, ρ_ref, 0.0, 5.0)
        x4_z8 = S.FourPosition(0.0, ρ_ref, 0.0, 8.0)
        x4_z2 = S.FourPosition(0.0, ρ_ref, 0.0, 2.0)

        @test p(geom, x4_z5) ≈ p(geom, x4_z8)
        @test p(geom, x4_z5) ≈ p(geom, x4_z2)
        @test p(geom, x4_z5) ≈ β_ref
    end
end

@testitem "RigidRotation profile with non-z-aligned axis" begin
    import Synchray as S
    using Synchray.Geometries

    # Conical geometry with axis along x-direction
    geom_x = Geometries.Conical(SVector(1, 0, 0), 0.1, 1.0 .. 10.0)

    β_ref = 0.2
    ρ_ref = 1.0
    p = Profiles.RigidRotation(β_ref, ρ_ref)

    @testset "ρ measured from x-axis" begin
        # Point at (5, 0, 1): distance from x-axis is 1
        x4_rho1 = S.FourPosition(0.0, 5.0, 0.0, 1.0)
        @test p(geom_x, x4_rho1) ≈ β_ref

        # Point at (5, 0, 2): distance from x-axis is 2
        x4_rho2 = S.FourPosition(0.0, 5.0, 0.0, 2.0)
        @test p(geom_x, x4_rho2) ≈ 2 * β_ref

        # Point at (5, 1, 1): distance from x-axis is √2
        x4_rho_sqrt2 = S.FourPosition(0.0, 5.0, 1.0, 1.0)
        @test p(geom_x, x4_rho_sqrt2) ≈ sqrt(2) * β_ref
    end

    @testset "On-axis gives zero" begin
        # Point on x-axis: (5, 0, 0)
        x4_on_axis = S.FourPosition(0.0, 5.0, 0.0, 0.0)
        @test p(geom_x, x4_on_axis) ≈ 0.0
    end
end

@testitem "CombinedVelocity" begin
    import Synchray as S
    using Synchray.Geometries, Synchray.Directions, Synchray.Profiles

    # Conical geometry
    geom = Geometries.Conical(SVector(0, 0, 1), 0.1, 1.0 .. 10.0)

    @testset "VelocitySpec addition creates CombinedVelocity" begin
        v1 = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(0.5))
        v2 = S.VelocitySpec(Directions.Toroidal(), S.beta, Profiles.Constant(0.3))

        combined = v1 + v2
        @test combined isa S.CombinedVelocity
        @test combined.v1 === v1
        @test combined.v2 === v2
    end

    @testset "Lab-frame β vector addition with hardcoded values" begin
        # Position (1, 0, 5) with z-axis cone:
        # - r = (1, 0, 5), |r| = √26
        # - e_r (radial from origin) = (1/√26, 0, 5/√26)
        # - e_φ (toroidal around z-axis) = (0, 1, 0)
        #
        # With β_r = 0.3, β_φ = 0.4:
        # - βv_r = 0.3 × (1/√26, 0, 5/√26)
        # - βv_φ = 0.4 × (0, 1, 0)
        # - βv_total = (0.3/√26, 0.4, 1.5/√26)
        # - |β_total|² = 0.09/26 + 0.16 + 2.25/26 = 0.09 + 0.16 = 0.25
        # - |β_total| = 0.5
        # - γ = 1/√(1-0.25) = 2/√3

        β_r = 0.3
        β_φ = 0.4
        v_radial = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(β_r))
        v_toroidal = S.VelocitySpec(Directions.Toroidal(), S.beta, Profiles.Constant(β_φ))
        combined = v_radial + v_toroidal

        x4 = S.FourPosition(0.0, 1.0, 0.0, 5.0)
        u_combined = S.four_velocity(combined, geom, x4)

        # Extract combined β vector
        βv_combined = @swiz(u_combined.xyz) / u_combined.t

        # Check individual components (hardcoded expected values)
        @test βv_combined[1] ≈ 0.3 / sqrt(26)  # x component
        @test βv_combined[2] ≈ 0.4             # y component
        @test βv_combined[3] ≈ 1.5 / sqrt(26)  # z component

        # Check total speed and Lorentz factor
        @test norm(βv_combined) ≈ 0.5
        @test u_combined.t ≈ 2 / sqrt(3)  # γ = 2/√3 ≈ 1.1547
    end

    @testset "CombinedVelocity with RigidRotation" begin
        # Radial + constant angular velocity
        β_r = 0.9
        β_φ_ref = 0.1
        ρ_ref = 1.0

        v_radial = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(β_r))
        v_rotation = S.VelocitySpec(Directions.Toroidal(), S.beta, Profiles.RigidRotation(β_φ_ref, ρ_ref))

        combined = v_radial + v_rotation

        # At ρ = ρ_ref: β_φ = β_φ_ref
        x4_at_ref = S.FourPosition(0.0, ρ_ref, 0.0, 5.0)
        u_at_ref = S.four_velocity(combined, geom, x4_at_ref)
        βv_at_ref = @swiz(u_at_ref.xyz) / u_at_ref.t
        @test norm(βv_at_ref) ≈ sqrt(β_r^2 + β_φ_ref^2) atol=1e-10

        # At ρ = 2 × ρ_ref: β_φ = 2 × β_φ_ref
        x4_at_2ref = S.FourPosition(0.0, 2*ρ_ref, 0.0, 5.0)
        u_at_2ref = S.four_velocity(combined, geom, x4_at_2ref)
        βv_at_2ref = @swiz(u_at_2ref.xyz) / u_at_2ref.t
        @test norm(βv_at_2ref) ≈ sqrt(β_r^2 + (2*β_φ_ref)^2) atol=1e-10
    end

    @testset "Full computation with unitful CombinedVelocity" begin
        using UnitfulAstro: pc

        # Create combined velocity: radial + unitful RigidRotation
        β_r = 0.3
        β_φ_ref = 0.4
        ρ_ref_unitful = 1.0pc

        v_radial = S.VelocitySpec(Directions.Radial(), S.beta, Profiles.Constant(β_r))
        v_rotation = S.VelocitySpec(Directions.Toroidal(), S.beta, Profiles.RigidRotation(β_φ_ref, ρ_ref_unitful))
        combined = v_radial + v_rotation

        # Strip units (as happens during actual computation)
        combined_stripped = S.ustrip(combined)

        # Create geometry (z-aligned cone)
        geom = S.Geometries.Conical(; axis=SVector(0, 0, 1), φj=0.1, z=1.0 .. 10.0)

        # Position with unitful coordinates, then convert to code units
        # At x = 1 pc from axis (ρ = ρ_ref), β_φ = β_φ_ref
        x4_unitful = S.FourPosition(0.0, 1.0, 0.0, 5.0) * pc
        x4_code = S._u_to_code(x4_unitful, S.UCTX.L0)

        # Compute velocity with unit-stripped spec and code-unit position
        u = S.four_velocity(combined_stripped, geom, x4_code)
        βv = @swiz(u.xyz) / u.t

        # Position (1, 0, 5) with z-axis cone:
        # - e_r (radial from origin) = (1/√26, 0, 5/√26)
        # - e_φ (toroidal around z-axis) = (0, 1, 0)
        # Expected: βv = 0.3 × e_r + 0.4 × e_φ = (0.3/√26, 0.4, 1.5/√26)
        @test βv[1] ≈ 0.3 / sqrt(26) atol=1e-10
        @test βv[2] ≈ 0.4 atol=1e-10
        @test βv[3] ≈ 1.5 / sqrt(26) atol=1e-10
        @test norm(βv) ≈ 0.5 atol=1e-10
    end
end

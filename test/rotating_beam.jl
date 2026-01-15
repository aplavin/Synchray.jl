@testitem "PrecessingNozzle - static case (θ_precession = 0)" begin
    import Synchray as S
    
    geom = S.Geometries.Conical(
        axis = SA[0.0, 0.0, 1.0],
        φj = 0.05,
        z = 1.0 .. 100.0
    )
    
    # With θ_precession = 0, nozzle axis always points along geometry axis
    # This should create a static central enhancement
    nozzle = S.Patterns.PrecessingNozzle(
        θ_precession = 0.0,
        θ_nozzle = 0.01,
        period = 100.0,
        β_flow = 1.0,
        profile = S.Patterns.TophatBump(5.0)
    )
    
    nozzle_prepared = S.prepare_for_computations(nozzle)
    
    # Test that unprepared and prepared nozzles give identical results
    x4_test = S.FourPosition(0.0, 0.0, 0.0, 10.0)
    @test nozzle(geom, x4_test, 1.0) ≈ nozzle_prepared(geom, x4_test, 1.0)
    
    # Test point on axis at different times - should always be enhanced
    for t in [0.0, 10.0, 50.0]
        x4 = S.FourPosition(t, 0.0, 0.0, 10.0)
        result = nozzle_prepared(geom, x4, 1.0)
        @test result ≈ 5.0
    end
    
    # Test point off-axis - should not be enhanced if outside nozzle cone
    x4_off = S.FourPosition(0.0, 1.0, 0.0, 10.0)
    result_off = nozzle_prepared(geom, x4_off, 1.0)
    @test result_off ≈ 1.0
end

@testitem "PrecessingNozzle - precessing case" begin
    import Synchray as S
    
    geom = S.Geometries.Conical(
        axis = SA[0.0, 0.0, 1.0],
        φj = 0.05,
        z = 1.0 .. 100.0
    )
    
    nozzle = S.Patterns.PrecessingNozzle(
        θ_precession = 0.02,  # 0.02 rad from axis
        θ_nozzle = 0.005,  # narrow nozzle
        period = 20.0,   # period
        φ0 = 0.0,
        β_flow = 1.0,
        profile = S.Patterns.TophatBump(10.0)
    )
    
    nozzle_prepared = S.prepare_for_computations(nozzle)
    
    # Specific hardcoded test: at t=0, r=10, check specific positions
    # Ejection time: t_ej = 0 - 10/1.0 = -10
    # Phase at t_ej = -10: 2π * (-10) / 20 = -π
    # At phase = -π: nozzle axis is at θ_precession from z-axis, pointing in -x direction
    # So nozzle axis ≈ [0, 0, cos(0.02)] + [sin(0.02) * cos(-π), 0, 0] ≈ [-0.02, 0, 0.9998]
    
    # Point directly in nozzle direction: should be maximally enhanced
    x4_in_beam = S.FourPosition(0.0, -0.2, 0.0, 10.0)
    result_in_beam = nozzle_prepared(geom, x4_in_beam, 1.0)
    @test result_in_beam ≈ 10.0
    
    # Point opposite to nozzle: should not be enhanced
    x4_opposite = S.FourPosition(0.0, 0.2, 0.0, 10.0)
    result_opposite = nozzle_prepared(geom, x4_opposite, 1.0)
    @test result_opposite ≈ 1.0
    
    # At t=5, r=10: t_ej = 5 - 10 = -5, phase = 2π*(-5)/20 = -π/2
    # Nozzle points in -y direction at phase = -π/2
    x4_t5_in_y = S.FourPosition(5.0, 0.0, -0.2, 10.0)
    result_t5 = nozzle_prepared(geom, x4_t5_in_y, 1.0)
    @test result_t5 ≈ 10.0
    
    # At t=0, nozzle should be in +x direction from axis
    # For a point at r=10, emission time is t_em = 0 - 10 = -10
    # Phase at t_em = -10: 2π * (-10) / 20 = -π
    # So nozzle points in -x direction at emission time
    
    # Point in -x direction should be enhanced
    x4_minus_x = S.FourPosition(0.0, -0.2, 0.0, 10.0)
    result = nozzle_prepared(geom, x4_minus_x, 1.0)
    # Check if this point is near the nozzle axis
    @test result > 1.0  # Should have some enhancement
    
    # Point far from nozzle axis should not be enhanced
    x4_far = S.FourPosition(0.0, 0.5, 0.5, 10.0)
    result_far = nozzle_prepared(geom, x4_far, 1.0)
    @test result_far ≈ 1.0
end

@testitem "PrecessingNozzle - time evolution" begin
    import Synchray as S
    
    geom = S.Geometries.Conical(
        axis = SA[0.0, 0.0, 1.0],
        φj = 0.05,
        z = 1.0 .. 100.0
    )
    
    nozzle = S.Patterns.PrecessingNozzle(
        θ_precession = 0.03,
        θ_nozzle = 0.01,
        period = 40.0,
        β_flow = 1.0,
        profile = S.Patterns.TophatBump(8.0)
    )
    
    nozzle_prepared = S.prepare_for_computations(nozzle)
    
    # Test that the same spatial point shows different enhancement
    # at different observation times due to rotation
    r = 20.0
    x_pos = SA[0.6, 0.0, r]
    
    results = Float64[]
    for t in [0.0, 10.0, 20.0, 30.0]
        x4 = S.FourPosition(t, x_pos...)
        push!(results, nozzle_prepared(geom, x4, 1.0))
    end
    
    # Results should vary due to rotation
    # (Not all should be the same)
    @test !all(r ≈ results[1] for r in results)
    
    # Specific hardcoded test: check specific spatial point at different times
    # Nozzle has θ_precession=0.03 rad, θ_nozzle=0.01 rad
    # At r=20: nozzle axis is at transverse distance ≈ 20 * sin(0.03) ≈ 0.6 from z-axis
    # Nozzle width at r=20: ≈ 20 * sin(0.01) ≈ 0.2
    # So points within 0.2 of the nozzle axis (at 0.6 transverse) should be enhanced
    
    # At t=0, r=20: t_ej = -20, phase = -π
    # Nozzle axis: [-sin(0.03), 0, cos(0.03)] ≈ [-0.03, 0, 1.0] (unnormalized direction)
    # At r≈20: nozzle center at x ≈ -0.6, nozzle width ≈ 0.2
    # Test point at x=-0.6 (nozzle center)
    x4_center_t0 = S.FourPosition(0.0, -0.6, 0.0, 20.0)
    result_center_t0 = nozzle_prepared(geom, x4_center_t0, 1.0)
    @test result_center_t0 ≈ 8.0
    
    # Same spatial point at t=10: phase = -π/2, nozzle in -y direction
    # Point at x=-0.6 is now far from nozzle
    result_center_t10 = nozzle_prepared(geom, S.FourPosition(10.0, -0.6, 0.0, 20.0), 1.0)
    @test result_center_t10 ≈ 1.0
    
    # At t=10: nozzle center should be at y ≈ -0.6
    result_y_t10 = nozzle_prepared(geom, S.FourPosition(10.0, 0.0, -0.6, 20.0), 1.0)
    @test result_y_t10 ≈ 8.0
end

@testitem "PrecessingNozzle - Integration with EmissionRegion" begin
    import Synchray as S
    
    # Create a full emission region with precessing nozzle pattern
    region = S.EmissionRegion(
        geometry = S.Geometries.Conical(
            axis = SA[0.0, 0.0, 1.0],
            φj = 0.05,
            z = 1.0 .. 50.0
        ),
        ne = S.Profiles.Modified(
            S.Profiles.Axial(S.PowerLaw(-2; val0=1e3, s0=1.0)),
            S.Patterns.PrecessingNozzle(
                θ_precession = 0.02,
                θ_nozzle = 0.005,
                period = 50.0,
                β_flow = 0.99,
                profile = S.Patterns.TophatBump(10.0)
            )
        ),
        B = S.BFieldSpec(
            S.Profiles.Axial(S.PowerLaw(-1; val0=1e-3, s0=1.0)),
            S.Directions.Scalar(),
            b -> S.FullyTangled(b)
        ),
        velocity = S.VelocitySpec(
            S.Directions.Axial(),
            S.Profiles.Constant(10.0)
        ),
        model = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=0.1)
    )
    
    region_prepared = S.prepare_for_computations(region)
    
    # Test that electron_density applies the pattern
    x4_test = S.FourPosition(0.0, 0.0, 0.0, 10.0)
    ne_val = S.electron_density(region_prepared, x4_test)
    
    # Should be a valid positive number
    @test ne_val > 0
    @test isfinite(ne_val)
end

@testitem "PrecessingNozzle - ustrip" begin
    import Synchray as S
    using Unitful
    
    nozzle = S.Patterns.PrecessingNozzle(
        θ_precession = 0.02,
        θ_nozzle = 0.005,
        period = 50.0u"s",
        β_flow = 0.99,
        profile = S.Patterns.TophatBump(10.0)
    )
    
    @test S.ustrip(nozzle) === S.Patterns.PrecessingNozzle(
        θ_precession = 0.02,
        θ_nozzle = 0.005,
        period = 1.49896229e12,
        β_flow = 0.99,
        profile = S.Patterns.TophatBump(10.0)
    )
end

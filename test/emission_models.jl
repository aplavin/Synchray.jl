@testitem "FixedEmission" begin
	import Synchray as S

	region = S.EmissionRegion(
		geometry = S.Geometries.Conical(; axis=SVector(0, 0, 1), φj=0.1, z=0.1..5),
		velocity = S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0)),
		emission = S.FixedEmission(S=1.0, α=0.01),
	)

	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=10.0, nz=256)

	I_val = S.render(ray, region)
	@test I_val > 0

	# frequency-independent: different k gives the same intensity
	ray2 = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=1.0, nz=256)
	@test S.render(ray2, region) ≈ I_val rtol=1e-10

	# higher S → higher intensity
	region_bright = @set region.emission.S = 10.0
	@test S.render(ray, region_bright) > I_val

	# higher α → more absorption → approaches S in thick limit
	region_thick = @set region.emission.α = 1e6
	@test S.render(ray, region_thick) ≈ 1.0 rtol=0.01
end

@testitem "PeakedEmission" begin
	import Synchray as S

	# use small α for optically thin regime (j = α · S · spectral is nonzero but thin)
	region = S.EmissionRegion(
		geometry = S.Geometries.Conical(; axis=SVector(0, 0, 1), φj=0.1, z=0.1..5),
		velocity = S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0)),
		emission = S.PeakedEmission(
			S = S.Profiles.Constant(1.0),
			α = S.Profiles.Constant(0.01),
			ν₀ = 10.0,
			σ = 0.3,
		),
	)

	ray_peak = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=10.0, nz=256)
	ray_off = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=1.0, nz=256)

	I_peak = S.render(ray_peak, region)
	@test I_peak > 0

	# off-peak intensity should be less
	@test S.render(ray_off, region) < I_peak

	# symmetry in log-frequency: ν₀*f and ν₀/f give equal intensity
	f = 2.0
	ray_hi = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=10.0 * f, nz=256)
	ray_lo = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=10.0 / f, nz=256)
	@test S.render(ray_hi, region) ≈ S.render(ray_lo, region) rtol=1e-10

	# optically thick limit: object should still be colored (peaked), not white (flat)
	# In the thick limit I_ν → S_ν = j_ν/α_ν, so S_ν must be spectrally peaked
	region_thick = @set region.emission.α = S.Profiles.Constant(1e6)
	I_thick_peak = S.render(ray_peak, region_thick)
	I_thick_off = S.render(ray_off, region_thick)
	@test I_thick_peak > 0
	@test I_thick_off < I_thick_peak / 2  # off-peak should be significantly dimmer

	# optically thick at peak: intensity saturates at source function value
	region_abs = @set region.emission.α = S.Profiles.Constant(10.0)
	@test S.render(ray_peak, region_abs) > I_peak  # thicker → closer to S_amp, higher than thin
	@test S.render(ray_peak, region_abs, S.OpticalDepth()) > 0

	# spectral index near peak should be close to zero
	si = S.render(ray_peak, region, S.SpectralIndex())
	@test abs(si) < 0.5

	# Float32 compatibility
	region32 = S.to_float_type(Float32, S.prepare_for_computations(region))
	ray32 = S.RayZ(; x0=S.FourPosition(0f0, 0f0, 0f0, 0f0), k=10f0, nz=64)
	@test S.render(ray32, region32) isa Float32
	@test S.render(ray32, region32) > 0
end

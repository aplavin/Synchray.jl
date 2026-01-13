@testitem "EmissionRegion geometry + scalings" begin
	import Synchray as S
	using Synchray: mean

	s0 = 1
	region = S.EmissionRegion(;
		geometry=S.Geometries.Conical(
			axis=SVector(0, 0, 1),
			φj=0.05,
			z=1e-3..5,
		),
		ne=S.Profiles.Axial(S.PowerLaw(-2; val0=2, s0)),
		B=S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=3, s0)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
		velocity=S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0)),
		model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1, Ca=1),
	)

	@testset "jet axis along ray" begin
		# Scalings along the axis at s=s0 and s=2s0.
		x4_1 = S.FourPosition(0, 0, 0, 1)
		x4_2 = S.FourPosition(0, 0, 0, 2)
		@test S.electron_density(region, x4_1) ≈ 2
		@test S.magnetic_field(region, x4_1) ≈ S.FullyTangled(3)
		@test S.electron_density(region, x4_2) ≈ 2 * (2^(-2))
		@test S.magnetic_field(region, x4_2) ≈ S.FullyTangled(3 * (2^(-1)))

		ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=1024)
		# On-axis ray should see the full truncation segment for axis=ẑ.
		@test S.z_interval(region, ray) == region.geometry.z
		@test S.render(ray, region) > 0
		# Off-axis ray within cone should also see a non-empty segment.
		xoff = 0.5 * maximum(region.geometry.z) * tan(region.geometry.φj)
		hit_ray = @set ray.x0.x = xoff
		@test S.z_interval(region, hit_ray) ≈ 2.5..5
		@test S.render(hit_ray, region) > 0

		# A ray outside the projected cone should miss entirely.
		xmiss = 1.1 * maximum(region.geometry.z) * tan(region.geometry.φj)
		miss_ray = @set ray.x0.x = xmiss
		@test S.z_interval(region, miss_ray) |> isempty
		@test S.render(miss_ray, region) == 0
	end
	@testset "off-axis viewing angles" begin
		s_probe = s0
		zmax = maximum(region.geometry.z)
		ymiss = 1.1 * zmax * tan(region.geometry.φj)

		cases = [
			(; label="small", θ=4 * region.geometry.φj, zpred = >(0)),
			(; label="large", θ=π / 4, zpred = >(0)),
			(; label="perpendicular", θ=π / 2, zpred = ∈(0 ± 1.1 * s_probe * tan(region.geometry.φj))),
			(; label="counterjet", θ=3π/4, zpred = <(0)),
		]

		@testset for (; label, θ, zpred) in cases
			regionθ = @set region.geometry.axis = SVector(sin(θ), 0, cos(θ))

			# Choose a ray that crosses the jet axis at s = s_probe.
			x0 = regionθ.geometry.axis.x * s_probe
			ray_zero = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=2048)
			ray_onaxis = @set ray_zero.x0.x = x0
			z_cross = regionθ.geometry.axis.z * s_probe
			x4_cross = S.FourPosition(0, x0, 0, z_cross)

			@test S.electron_density(regionθ, x4_cross) > 0

			zint = S.z_interval(regionθ, ray_zero)
			@test isempty(zint)
			zint = S.z_interval(regionθ, ray_onaxis)
			@test !isempty(zint)
			@test all(zpred, endpoints(zint))
			@test z_cross ∈ zint
			@test S.render(ray_onaxis, regionθ) > 0

			ray_half = @set ray_onaxis.x0.y = 0.5 * s_probe * tan(region.geometry.φj)
			zint_half = S.z_interval(regionθ, ray_half)
			@test mean(zint_half) ≈ mean(zint)
			@test width(zint_half) ≈ √(1 - 0.5^2) * width(zint)  rtol=2e-2

			# A ray with |y| larger than the maximal jet radius must miss, regardless of tilt.
			ray_miss = @set ray_onaxis.x0.y = ymiss
			@test S.z_interval(regionθ, ray_miss) |> isempty
			@test S.render(ray_miss, regionθ) == 0
			@test S.render(ray_miss, regionθ, S.OpticalDepth()) == 0
			@test S.render(ray_miss, regionθ, S.SpectralIndex()) |> isnan

			ray_miss = @modify(-, ray_onaxis.x0.x)
			@test S.z_interval(regionθ, ray_miss) |> isempty
			@test S.render(ray_miss, regionθ) == 0
		end
	end
end


# @testitem "Unitful boundary API (EmissionRegion)" begin
# 	import Synchray as S
# 	using Unitful, UnitfulAstro
# 	using RectiGrids

# 	region = S.withunits(S.EmissionRegion;
# 		geometry=S.Geometries.Conical(
# 			axis=SVector(0, 0, 1),
# 			φj=2u"°",
# 			z=1e-3u"pc"..50u"pc",
# 		),
# 		ne=S.Profiles.Axial(S.PowerLaw(-2; val0=2u"cm^-3", s0=1u"pc")),
# 		B=S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=3u"Gauss", s0=1u"pc")), S.Directions.Scalar(), b -> S.FullyTangled(b)),
# 		velocity=S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0.0)),
# 		model=S.IsotropicPowerLawElectrons(; p=2.5),
# 	)

# 	@test region.model.Cj_ordered * region.model.sinavg_j ≈ 1.7e-18 rtol=1e-3
# 	@test region.model.Ca_ordered * region.model.sinavg_a ≈ 7.20e13 rtol=1e-3

# 	cam = S.withunits(S.CameraZ;
# 		xys=grid(SVector, range(0.01u"pc"..0.1u"pc", 2), range(-0.001u"pc"..0.001u"pc", 2)),
# 		nz=20,
# 		ν=230u"GHz",
# 		t=0u"yr",
# 	)

# 	Iν_nou = S.render(cam, region)
# 	@test all(>(0), Iν_nou)
# 	Iν = S.withunits(S.render, cam, region, S.Intensity())
# 	@test unit(eltype(Iν)) == u"erg/s/cm^2/Hz/sr"
# 	@test ustrip.(Iν) == Iν_nou
# 	@test ustrip.(u"Jy/mas^2", Iν) ≈ [121 121; 0.028 0.028]  rtol=1e-2

# 	res = S.withunits(S.render, cam, region, (S.Intensity(), S.SpectralIndex()))
# 	@test res[1] == Iν
# 	@test eltype(res[2]) == Float64
# 	@test size(res[2]) == size(Iν)
# end


@testitem "EmissionRegion phenomenology" begin
	import Synchray as S
	using Accessors
	using RectiGrids
	using Optim, Roots

	φj = 0.05
	θ = 0.2  # viewing angle (> φj)
	s0 = 1.

	region = S.EmissionRegion(;
		geometry=S.Geometries.Conical(;
			axis=S.SVector(sin(θ), 0.0, cos(θ)),
			φj,
			z=1e-3..50,
		),
		ne=S.Profiles.Axial(S.PowerLaw(-2; val0=1., s0)),
		B=S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=3., s0)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
		velocity=S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0f0)),
		model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1.0),
	)

	@testset "core shift (intensity peak along axis) scales ~ ν^-1" begin
		core_byint(region, ν) = begin
			ray_base = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=ν, nz=2048)
			opt = Optim.optimize(x -> -S.render((@set ray_base.x0.x = x), region), 0.001, 5)
			return Optim.minimizer(opt)
		end
		core_byτ(region, ν) = begin
			ray_base = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=ν, nz=2048)
			Roots.find_zero(x -> S.render((@set ray_base.x0.x = x), region, S.OpticalDepth()) - 1, (0.001, 5))
		end

		@test core_byint(region, 1) ≈ 2*core_byint(region, 2) rtol=1e-2
		@test core_byint(region, 4) ≈ 0.5*core_byint(region, 2) rtol=1e-2
		@test core_byτ(region, 1) ≈ 2*core_byτ(region, 2) rtol=1e-2
		@test core_byτ(region, 4) ≈ 0.5*core_byτ(region, 2) rtol=1e-2
		@test 0.7core_byτ(region, 1) ≈ core_byint(region, 1) rtol=0.05
	end

	@testset "spectrum: integrated flux density is ~flat" begin
		# Integrate over a window that contains the full truncated jet projection.
		cam0 = S.CameraZ(; xys=grid(S.SVector, range(-0.1..20, 201), range(-5..5, 201)), nz=20, ν=NaN, t=0)
		flux(ν) = begin
			cam = @set cam0.ν = ν
			img = S.render(cam, region)
			@assert iszero(img[end,:])
			@assert iszero(img[:,end])
			sum(img)
		end

		@test flux(1) > 0
		@test flux(1) / flux(2) ≈ 1  rtol=0.02
		@test flux(0.5) / flux(1) ≈ 1  rtol=0.07
	end

	ray_at_s(ν, s) = begin
		rxy = (@swiz region.geometry.axis.xy) * s
		S.RayZ(; x0=S.FourPosition(0, rxy..., 0), k=ν, nz=4096)
	end

	@testset "thin-regime scaling with (ne0, B0) matches IsotropicPowerLawElectrons" begin
		νthin = 80.0
		ray = ray_at_s(νthin, s0)
		τthin = S.render(ray, region, S.OpticalDepth())
		@test τthin < 0.2

		I0 = S.render(ray, region)
		I_ne = S.render(ray, @set region.ne.f.val0 *= 2)
		@test I_ne ≈ 2I0 rtol=0.05

		p = region.model.p
		I_B = S.render(ray, @set region.B.scale.f.val0 *= 2)
		@test I_B ≈ I0 * (2^((p + 1) / 2)) rtol=0.07
	end
end

@testitem "float64 vs float32" begin
	import Synchray as S
	using Accessors
	using RectiGrids

	φj = 0.05
	θ = 0.2  # viewing angle (> φj)

	region = S.EmissionRegion(;
		geometry=S.Geometries.Conical(;
			axis=S.SVector(sin(θ), 0.0, cos(θ)),
			φj,
			z=1e-3..50,
		),
		ne=S.Profiles.Axial(S.PowerLaw(-2; val0=1., s0=1.)),
		B=S.BFieldSpec(S.Profiles.Axial(S.PowerLaw(-1; val0=3., s0=1.)), S.Directions.Scalar(), b -> S.FullyTangled(b)),
		velocity=S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0f0)),
		model=S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1.0),
	)

	cam64 = S.CameraZ(; xys=grid(S.SVector, range(0.01..0.1, 2), range(-0.001..0.001, 2)), nz=20, ν=2., t=0.)
	cam32 = S.to_float_type(Float32, cam64)

	region64 = region
	region32 = S.to_float_type(Float32, S.prepare_for_computations(region64))

	f64 = S.render(cam64, region64)
	f32 = S.render(cam32, region32)
	@test eltype(f64) == Float64
	@test eltype(f32) == Float32
	@test all(>(0), f64)
	@test f32 ≈ f64  rtol=1e-5

	# Test that coordinate transformation functions preserve Float32
	rot_mat32 = S.rotation_lab_to_local(region32.geometry)
	rot_mat64 = S.rotation_lab_to_local(region64.geometry)
	@test eltype(rot_mat32) == Float32
	@test eltype(rot_mat64) == Float64
	@test rot_mat32 ≈ rot_mat64

	r32 = S.SVector(Float32(1.5), Float32(-0.3), Float32(2.0))
	r64 = S.SVector(1.5, -0.3, 2.0)

	rlocal32 = S.rotate_lab_to_local(region32.geometry, r32)
	rlocal64 = S.rotate_lab_to_local(region64.geometry, r64)
	@test eltype(rlocal32) == Float32
	@test eltype(rlocal64) == Float64
	@test rlocal32 ≈ rlocal64

	rlab32 = S.rotate_local_to_lab(region32.geometry, rlocal32)
	rlab64 = S.rotate_local_to_lab(region64.geometry, rlocal64)
	@test eltype(rlab32) == Float32
	@test eltype(rlab64) == Float64
	@test rlab32 ≈ r32
	@test rlab64 ≈ r64
end

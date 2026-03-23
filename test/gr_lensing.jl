@testitem "GR lensing" begin
	import Synchray as S
	using Krang

	spin = 0.5
	θ_obs = π / 3  # 60° inclination

	cam = S.CameraOrtho(;
		look_direction = S.SVector(0.0, 0.0, -1.0),
		xys = S.grid(S.SVector, x=range(-60.0, 60.0, length=32), y=range(-60.0, 60.0, length=32)),
		nz = 100, ν = 1e11, t = 0.0,
	)

	defl = S.compute_deflection_map(spin, θ_obs, cam)

	# Incoming asymptotic direction in BL→Cartesian coords (observer at θ_obs from spin axis)
	n_in = S.SVector(0.0, sin(θ_obs), -cos(θ_obs))

	@testset "shadow detection" begin
		for I in CartesianIndices(defl)
			b = S.norm(cam.xys[I])
			if b < 3.0
				@test defl[I] === nothing
			end
			if b > 15.0
				@test defl[I] !== nothing
			end
		end
	end

	@testset "outgoing direction is unit vector" begin
		for d in defl
			isnothing(d) && continue
			@test S.norm(S.direction3(d)) ≈ 1 atol=1e-10
		end
	end

	@testset "anchor at large r" begin
		for d in defl
			isnothing(d) && continue
			@test S.norm(S.SVector(d.x0.x, d.x0.y, d.x0.z)) > 10.0
		end
	end

	@testset "weak-field deflection angle" begin
		# For b >> M, deflection ≈ 4M/b (Schwarzschild limit, M=1 in Krang)
		# Deflection = angle between incoming and outgoing asymptotic directions
		for I in CartesianIndices(defl)
			d = defl[I]
			isnothing(d) && continue
			b = S.norm(cam.xys[I])
			b < 40.0 && continue

			# For zero deflection, n_out = n_in (straight through).
			# Deflection angle = acos(dot(n_in, n_out))
			defl_angle = acos(clamp(S.dot(n_in, S.direction3(d)), -1, 1))
			predicted = 4.0 / b
			@test defl_angle ≈ predicted rtol=0.25
		end
	end

	@testset "scattering plane" begin
		# For β≈0 pixels, the scattering plane contains the incoming direction and the
		# x-axis offset. Its normal is cross(x̂, n_in) ∝ (0, cos θ, sin θ).
		expected_normal = S.normalize(S.SVector(0.0, cos(θ_obs), sin(θ_obs)))
		for I in CartesianIndices(defl)
			d = defl[I]
			isnothing(d) && continue
			uv = cam.xys[I]
			abs(uv[2]) > 0.1 && continue
			abs(uv[1]) < 15.0 && continue

			plane_normal = S.normalize(S.cross(n_in, S.direction3(d)))
			@test abs(S.dot(plane_normal, expected_normal)) > 0.9
		end
	end

	@testset "GR render ≥ flat render" begin
		region = S.EmissionRegion(;
			geometry = S.Geometries.Conical(; axis=S.SVector(0.0, 0.0, 1.0), φj=deg2rad(15), z=20.0..100.0),
			velocity = S.VelocitySpec(S.Directions.Axial(), S.gamma, S.Profiles.Constant(2.0)),
			emission = S.SynchrotronEmission(
				ne = S.Profiles.Axial(S.Profiles.PowerLaw(-1.5; val0=100.0, s0=1.0)),
				B = S.BFieldSpec(
					S.Profiles.Axial(S.Profiles.PowerLaw(-1.8; val0=500.0, s0=1.0)),
					S.Directions.Scalar(),
					b -> S.FullyTangled(b)
				),
				electrons = S.IsotropicPowerLawElectrons(; p=3.0, γmin=100.0, γmax=1e5),
			),
		) |> S.prepare_for_computations

		cam_small = S.CameraOrtho(;
			look_direction = S.SVector(0.0, 0.0, -1.0),
			xys = S.grid(S.SVector, x=range(-30.0, 30.0, length=8), y=range(-30.0, 30.0, length=8)),
			nz = 100, ν = 1e11, t = 0.0,
		)
		defl_small = S.compute_deflection_map(spin, θ_obs, cam_small)

		img_flat = S.render(cam_small, region)
		img_gr = S.render(S.CameraGR(; camera=cam_small, deflection=defl_small), region)

		@test sum(img_gr) ≥ sum(img_flat) - 1e-30
		@test img_gr != img_flat
	end
end

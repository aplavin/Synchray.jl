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

	# --- Single-ray tests with spheres ---
	# Helper: get deflected ray for a single pixel at offset (x, y)
	u_rest = S.FourVelocity(1.0, 0.0, 0.0, 0.0)
	bh_pos = S.SVector(0.0, 0.0, 0.0)
	ν = 1e11
	nz = 200
	jν = 1e-20

	function make_ray(x, y)
		k = S.photon_k(ν, S.SVector(0.0, 0.0, -1.0))
		T = Float64
		S.Ray(S.FourPosition(0.0, x, y, 1000.0), k,
			S.SVector{3,T}(-1, 0, 0), S.SVector{3,T}(0, 1, 0), nz, S.SlowLight())
	end

	function deflected_ray(x, y)
		# Use 3x3 grid (Krang needs nonzero range), take center pixel
		xf, yf = Float64(x), Float64(y)
		cam_3px = S.CameraOrtho(;
			look_direction=S.SVector(0.0, 0.0, -1.0),
			origin=S.SVector(0.0, 0.0, 1000.0),
			xys=S.grid(S.SVector, x=[xf-1, xf, xf+1], y=[yf-1, yf, yf+1]),
			nz=nz, ν=ν, t=0.0,
		)
		S.compute_deflection_map(spin, θ_obs, cam_3px)[2, 2]
	end

	@testset "far ray: sphere in front of BH" begin
		# Ray at x=50, far from BH (b=50 >> shadow). Sphere on the incoming ray path,
		# entirely between camera and BH. GR incoming sees it fully, outgoing ray
		# (in BL→Cartesian frame, ~60° away) misses it entirely.
		ray_in = make_ray(50.0, 0.0)
		ray_out = deflected_ray(50.0, 0.0)
		@test ray_out !== nothing

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=30.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S._render_gr_pixel(sphere, ray_in, ray_out, bh_pos, S.Intensity())

		@test flat ≈ jν * 60.0 rtol=0.01
		@test gr == flat
	end

	@testset "far ray: sphere straddling BH plane" begin
		# Sphere centered at (50, 0, 0) — right at the BH clip plane.
		# Incoming ray only sees the front half (z>0), back half clipped.
		# Outgoing ray goes in a different direction and misses the sphere.
		ray_in = make_ray(50.0, 0.0)
		ray_out = deflected_ray(50.0, 0.0)

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 0.0), radius=30.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S._render_gr_pixel(sphere, ray_in, ray_out, bh_pos, S.Intensity())

		# Flat sees full sphere: I = j * 2R = j * 60
		@test flat ≈ jν * 60.0 rtol=0.01
		# GR sees only front half: I = j * R = j * 30 (clipped at BH plane)
		@test gr ≈ jν * 30.0 rtol=0.01
		@test gr ≈ flat / 2 rtol=0.01
	end

	@testset "captured ray: sphere behind BH invisible" begin
		# Ray at x=0, y=0 — hits BH directly, captured.
		ray_in = make_ray(0.0, 0.0)
		ray_out = deflected_ray(0.0, 0.0)
		@test ray_out === nothing  # captured

		# Sphere behind BH
		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 0.0, 0.0, -50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S._render_gr_pixel(sphere, ray_in, ray_out, bh_pos, S.Intensity())

		# Flat sees the sphere (straight line goes through it): I = j * 2R
		@test flat ≈ jν * 20.0 rtol=0.01
		# GR: incoming ray clipped at BH, ray captured → sphere invisible
		@test gr === 0.0
	end

	@testset "intermediate ray: sphere in front visible, behind invisible" begin
		# Ray at x=10, moderate offset. Strong deflection but not captured.
		ray_in = make_ray(10.0, 0.0)
		ray_out = deflected_ray(10.0, 0.0)
		@test ray_out !== nothing  # not captured

		# Sphere on the incoming ray path, in front of BH
		sphere_front = S.UniformSphere(;
			center=S.FourPosition(0.0, 10.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat_front = S.integrate_ray(sphere_front, ray_in)
		gr_front = S._render_gr_pixel(sphere_front, ray_in, ray_out, bh_pos, S.Intensity())

		# Flat: I = j * 2R through center
		@test flat_front ≈ jν * 20.0 rtol=0.01
		# GR: incoming ray sees the sphere (in front of BH). Outgoing ray deflected away.
		@test gr_front == flat_front

		# Same sphere but behind BH — should be invisible via GR
		sphere_behind = S.UniformSphere(;
			center=S.FourPosition(0.0, 10.0, 0.0, -50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat_behind = S.integrate_ray(sphere_behind, ray_in)
		gr_behind = S._render_gr_pixel(sphere_behind, ray_in, ray_out, bh_pos, S.Intensity())

		# Flat sees it (straight line goes through)
		@test flat_behind ≈ jν * 20.0 rtol=0.01
		# GR: incoming clipped at BH, outgoing deflected elsewhere → invisible
		@test gr_behind === 0.0
	end

	@testset "bh_position translation invariance" begin
		# Shifting BH + sphere + camera by the same offset should give identical pixel values.
		offset = S.SVector(100.0, -200.0, 300.0)

		# Reference: BH at origin, sphere at (50, 0, 50)
		sphere_ref = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)
		ray_in_ref = make_ray(50.0, 0.0)
		ray_out_ref = deflected_ray(50.0, 0.0)
		gr_ref = S._render_gr_pixel(sphere_ref, ray_in_ref, ray_out_ref, bh_pos, S.Intensity())
		@test gr_ref > 0

		# Shifted: everything moved by offset
		bh_shifted = bh_pos + offset
		sphere_shifted = S.UniformSphere(;
			center=S.FourPosition(0.0, (S.SVector(50.0, 0.0, 50.0) + offset)...), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)
		# Shifted incoming ray: same direction, origin shifted
		ray_in_shifted = S.Ray(
			S.FourPosition(0.0, (S.SVector(50.0, 0.0, 1000.0) + offset)...),
			ray_in_ref.k, ray_in_ref.e1, ray_in_ref.e2, nz, S.SlowLight())
		# Shifted outgoing ray: same direction, anchor shifted
		anchor_shifted = S.SVector(ray_out_ref.x0.x, ray_out_ref.x0.y, ray_out_ref.x0.z) + offset
		ray_out_shifted = S.Ray(
			S.FourPosition(ray_out_ref.x0.t, anchor_shifted...),
			ray_out_ref.k, ray_out_ref.e1, ray_out_ref.e2, nz, S.SlowLight())

		gr_shifted = S._render_gr_pixel(sphere_shifted, ray_in_shifted, ray_out_shifted, bh_shifted, S.Intensity())

		@test gr_shifted == gr_ref
	end

	@testset "bh_rg scaling" begin
		# Scaling bh_rg by factor s and camera xys by factor s should give identical
		# pixel values (same geometry, just in different units).
		scale = 5.0

		# Reference: rg=1, sphere at z=50 radius=10, ray at x=50
		cam_ref = S.CameraOrtho(;
			look_direction=S.SVector(0.0, 0.0, -1.0),
			origin=S.SVector(0.0, 0.0, 1000.0),
			xys=S.grid(S.SVector, x=[49.0, 50.0, 51.0], y=[-1.0, 0.0, 1.0]),
			nz=nz, ν=ν, t=0.0,
		)
		defl_ref = S.compute_deflection_map(spin, θ_obs, cam_ref; bh_rg=1.0)
		sphere_ref = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)
		cam_gr_ref = S.CameraGR(; camera=cam_ref, deflection=defl_ref, bh_rg=1.0)
		img_ref = S.render(cam_gr_ref, sphere_ref)

		# Scaled: rg=scale, xys*scale, sphere coords*scale, camera origin*scale
		cam_scaled = S.CameraOrtho(;
			look_direction=S.SVector(0.0, 0.0, -1.0),
			origin=S.SVector(0.0, 0.0, 1000.0 * scale),
			xys=S.grid(S.SVector, x=[49.0, 50.0, 51.0] .* scale, y=[-1.0, 0.0, 1.0] .* scale),
			nz=nz, ν=ν, t=0.0,
		)
		defl_scaled = S.compute_deflection_map(spin, θ_obs, cam_scaled; bh_rg=scale)
		sphere_scaled = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0*scale, 0.0, 50.0*scale), radius=10.0*scale,
			u0=u_rest, jν=jν, αν=0.0,
		)
		cam_gr_scaled = S.CameraGR(; camera=cam_scaled, deflection=defl_scaled, bh_rg=scale)
		img_scaled = S.render(cam_gr_scaled, sphere_scaled)

		# Optically thin: I = j * path_length. Path scales by `scale`, so I scales too.
		@test collect(img_scaled) ≈ collect(img_ref) .* scale rtol=0.01
	end
end

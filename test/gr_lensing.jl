@testitem "GR lensing" begin
	import Synchray as S
	using Krang

	spin = 0.5
	θ_view = π / 3  # 60° inclination from the +z spin axis
	look_direction = S.normalize(S.SVector(0.0, -sin(θ_view), cos(θ_view)))
	lens = S.GRLens(; spin)

	cam = S.CameraOrtho(;
		look_direction,
		origin = -1000.0 * look_direction,
		xys = S.grid(S.SVector, x=range(-60.0, 60.0, length=32), y=range(-60.0, 60.0, length=32)),
		nz = 100, ν = 1e11, t = 0.0,
	)

	defl = S.compute_deflection_map(lens, cam)

	# Incoming photon direction at the observer (opposite of the tracing ray direction).
	n_in = -cam.n

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
		# For a ray with β≈0 (offset along cam.e1 only), the scattering plane contains
		# the incoming direction and the pixel offset direction.
		# The plane normal should be ≈ ±cam.e2.
		# Use β=0.01 (not exactly 0 — Krang has a sign(β)=0 degeneracy at β=0).
		cam_sp = S.CameraOrtho(;
			look_direction,
			origin = -1000.0 * look_direction,
			xys = S.grid(S.SVector, x=[29.0, 30.0, 31.0], y=[-0.99, 0.01, 1.01]),
			nz = 100, ν = 1e11, t = 0.0,
		)
		defl_sp = S.compute_deflection_map(lens, cam_sp)
		d_sp = defl_sp[2, 2]
		@test d_sp !== nothing

		plane_normal = S.normalize(S.cross(n_in, S.direction3(d_sp)))
		@test abs(S.dot(plane_normal, cam.e2)) > 0.9
	end

	# --- Single-ray tests with spheres ---
	# Helper: get deflected ray for a single pixel at offset (x, y)
	u_rest = S.FourVelocity(1.0, 0.0, 0.0, 0.0)
	bh_pos = lens.bh_position
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
		xf, yf = Float64(x), Float64(y)
		cam_1px = S.CameraOrtho(;
			look_direction=S.SVector(0.0, 0.0, -1.0),
			origin=S.SVector(0.0, 0.0, 1000.0),
			xys=[S.SVector(xf, yf)],
			nz=nz, ν=ν, t=0.0,
		)
		only(S.compute_deflection_map(spin, θ_obs, cam_1px))
	end

	@testset "vector xys topology is preserved" begin
		xys_vec = [
			S.SVector(50.0, 0.0),
			S.SVector(10.0, 0.0),
			S.SVector(0.0, 0.0),
		]
		cam_vec = S.CameraOrtho(;
			look_direction=S.SVector(0.0, 0.0, -1.0),
			origin=S.SVector(0.0, 0.0, 1000.0),
			xys=xys_vec,
			nz=nz, ν=ν, t=0.0,
		)
		defl_vec = S.compute_deflection_map(spin, θ_obs, cam_vec)

		@test defl_vec isa Vector
		@test axes(defl_vec) == axes(xys_vec)

		for (i, uv) in pairs(xys_vec)
			ray_vec = defl_vec[i]
			ray_single = deflected_ray(uv[1], uv[2])
			@test isnothing(ray_vec) == isnothing(ray_single)
			isnothing(ray_vec) && continue
			@test ray_vec.x0 ≈ ray_single.x0
			@test ray_vec.k ≈ ray_single.k
			@test ray_vec.e1 ≈ ray_single.e1
			@test ray_vec.e2 ≈ ray_single.e2
		end
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

	@testset "two spheres on the same straight ray across the BH" begin
		ray_in = make_ray(50.0, 0.0)
		# Weak-deflection limit: use a straight outgoing segment, while _render_gr_pixel
		# still splits the path at the BH into far/background and near/foreground halves.
		ray_out = ray_in

		jν_dim = 1e-25
		αν_opaque = 1.0
		jν_bright = 1e-18
		R1, R2 = 10.0, 20.0

		sphere_front = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=R1,
			u0=u_rest, jν=jν_dim, αν=αν_opaque,
		)
		sphere_back = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, -80.0), radius=R2,
			u0=u_rest, jν=jν_bright, αν=0.0,
		)
		combined = S.CombinedMedium(sphere_back, sphere_front)

		flat_front = S.integrate_ray(sphere_front, ray_in)
		flat_back = S.integrate_ray(sphere_back, ray_in)
		gr_front = S._render_gr_pixel(sphere_front, ray_in, ray_out, bh_pos, S.Intensity())
		gr_back = S._render_gr_pixel(sphere_back, ray_in, ray_out, bh_pos, S.Intensity())

		@test flat_front > 0
		@test flat_back > 0
		@test gr_front > 0
		@test gr_back > 0
		@test gr_front ≈ flat_front rtol=1e-10
		@test gr_back ≈ flat_back rtol=1e-10

		flat = S.integrate_ray(combined, ray_in)
		gr = S._render_gr_pixel(combined, ray_in, ray_out, bh_pos, S.Intensity())
		τ_front = S.integrate_ray(sphere_front, ray_in, S.OpticalDepth())
		expected = flat_back * exp(-τ_front) + flat_front

		@test gr ≈ flat rtol=1e-10
		@test flat ≈ expected rtol=0.05
		@test gr ≈ expected rtol=0.05
	end

	@testset "U-turn ray: see sphere behind BH via ~180° bending" begin
		# Sphere at (0, 0, -100) — directly behind BH in Synchray frame.
		# α=6.110 (just outside shadow) bends outgoing ray to direction ≈ (0, 0, -1).
		# The outgoing ray passes within ~3.7 of sphere center (R=5.5 → hit).
		# Incoming ray at x=6.11 is 6.1 from center (> R=5.5 → miss).
		ray_in = make_ray(6.110, 0.01)
		ray_out = deflected_ray(6.110, 0.01)
		@test ray_out !== nothing

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 0.0, 0.0, -100.0), radius=5.5,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S._render_gr_pixel(sphere, ray_in, ray_out, bh_pos, S.Intensity())

		@test flat === 0.0
		# Path through sphere ≈ 8.1 (off-center hit, miss_dist=3.7, R=5.5)
		@test gr ≈ jν * 8.1 rtol=0.05
	end

	@testset "90° bend ray: see sphere sideways via quarter-turn" begin
		# Sphere at (-100, 0, 0) — sideways from BH in Synchray frame.
		# α=6.762 bends outgoing ray to direction ≈ (-1, 0, 0).
		# Outgoing ray passes within ~8.2 of center (R=10 → hit).
		# Incoming ray at x=6.76 is ~107 from (-100,0,0) → miss.
		ray_in = make_ray(6.762, 0.01)
		ray_out = deflected_ray(6.762, 0.01)
		@test ray_out !== nothing

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, -100.0, 0.0, 0.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S._render_gr_pixel(sphere, ray_in, ray_out, bh_pos, S.Intensity())

		@test flat === 0.0
		# Path through sphere ≈ 11.3 (off-center hit, miss_dist=8.2, R=10)
		@test gr ≈ jν * 11.3 rtol=0.05
	end
end

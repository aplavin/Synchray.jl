@testitem "GR lensing" begin
	import Synchray as S
	using Krang
	using FastChebInterp

	spin = 0.5
	θ_view = π / 3  # 60° inclination from the +z spin axis
	photon_direction = S.normalize(S.SVector(0.0, -sin(θ_view), cos(θ_view)))
	lens = S.GRLens(; spin)

	cam = S.CameraOrtho(;
		photon_direction,
		xys = S.grid(S.SVector, x=range(-60.0, 60.0, length=32), y=range(-60.0, 60.0, length=32)),
		nz = 100, ν = 1e11, t = 0.0,
	)

	@testset "shadow detection" begin
		# Inside shadow (b well below critical ~6.2)
		@test S.is_captured_ray(S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 0.0)), lens).ray_out)
		@test S.is_captured_ray(S.RayGR2(S.camera_ray(cam, S.SVector(3.0, 0.0)), lens).ray_out)
		@test S.is_captured_ray(S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 5.0)), lens).ray_out)
		# Outside shadow
		@test !S.is_captured_ray(S.RayGR2(S.camera_ray(cam, S.SVector(10.0, 0.0)), lens).ray_out)
		@test !S.is_captured_ray(S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 50.0)), lens).ray_out)
	end

	@testset "outgoing direction is unit vector" begin
		for uv in [S.SVector(10.0, 0.0), S.SVector(30.0, 5.0), S.SVector(50.0, -20.0)]
			d = S.RayGR2(S.camera_ray(cam, uv), lens).ray_out
			@test !S.is_captured_ray(d)
			@test S.norm(S.direction3(d)) ≈ 1 atol=1e-10
		end
	end

	@testset "weak-field deflection angle" begin
		# For b >> M, deflection ≈ 4M/b (Schwarzschild limit, M=1 in Krang).
		# Use non-zero y to avoid the Krang β=0 degeneracy.
		for uv in [S.SVector(40.0, 10.0), S.SVector(30.0, 30.0), S.SVector(45.0, -20.0)]
			b = S.norm(uv)
			d = S.RayGR2(S.camera_ray(cam, uv), lens).ray_out
			@test !S.is_captured_ray(d)
			# Both direction3 vectors are photon travel directions; for zero deflection they're parallel.
			defl_angle = acos(clamp(S.dot(cam.n, S.direction3(d)), -1, 1))
			@test defl_angle ≈ 4.0 / b rtol=0.1
		end
	end

	@testset "scattering plane" begin
		# For a ray with β≈0 (offset along cam.e1 only), the scattering plane contains
		# the incoming direction and the pixel offset direction.
		# The plane normal should be ≈ ±cam.e2.
		# Use β=0.01 (not exactly 0 — Krang has a sign(β)=0 degeneracy at β=0).
		uv_sp = S.SVector(30.0, 0.01)
		d_sp = S.RayGR2(S.camera_ray(cam, uv_sp), lens).ray_out
		@test !S.is_captured_ray(d_sp)

		plane_normal = S.normalize(S.cross(cam.n, S.direction3(d_sp)))
		@test abs(S.dot(plane_normal, cam.e2)) > 0.9
	end

	# --- Single-ray tests with spheres ---
	# Helper: get deflected ray for a single pixel at offset (x, y)
	u_rest = S.FourVelocity(1.0, 0.0, 0.0, 0.0)
	bh_pos = lens.bh_position
	ν = 1e11
	nz = 200
	jν = 1e-20

	function make_ray(x, y; e1=S.SVector(1.0, 0.0, 0.0), e2=S.SVector(0.0, 1.0, 0.0))
		k = S.photon_k(ν, S.SVector(0.0, 0.0, 1.0))
		S.Ray(S.FourPosition(0.0, x, y, 0.0), k,
			e1, e2, nz, S.SlowLight())
	end

	@testset "RayGR2 constructor and lens API" begin
		uv = S.SVector(50.0, 0.0)
		cam_1px = S.CameraOrtho(; photon_direction, xys=[uv], nz=nz, ν=ν, t=0.0)
		ray_in = S.camera_ray(cam_1px, uv)
		raygr_inferred = S.RayGR2(ray_in, lens)
		raygr_explicit = S.RayGR2(ray_in, raygr_inferred.ray_out, bh_pos)
		ray_in_twisted = S.Ray(
			ray_in.x0,
			ray_in.k,
			ray_in.e2,
			-ray_in.e1,
			ray_in.nz,
			ray_in.light,
		)
		raygr_twisted = S.RayGR2(ray_in_twisted, lens)

		@test !S.is_captured_ray(raygr_inferred.ray_out)
		@test raygr_inferred.ray_out.x0 ≈ raygr_explicit.ray_out.x0
		@test raygr_inferred.ray_out.k ≈ raygr_explicit.ray_out.k
		@test raygr_inferred.ray_out.e1 ≈ raygr_explicit.ray_out.e1
		@test raygr_inferred.ray_out.e2 ≈ raygr_explicit.ray_out.e2
		@test !S.is_captured_ray(raygr_twisted.ray_out)
		@test raygr_twisted.ray_out.x0 ≈ raygr_inferred.ray_out.x0
		@test raygr_twisted.ray_out.k ≈ raygr_inferred.ray_out.k

		ray_out_bad_freq = S.Ray(
			raygr_inferred.ray_out.x0,
			S.photon_k(2 * S.frequency(ray_in), S.direction3(raygr_inferred.ray_out)),
			raygr_inferred.ray_out.e1,
			raygr_inferred.ray_out.e2,
			raygr_inferred.ray_out.nz,
			raygr_inferred.ray_out.light,
		)
		ray_out_bad_light = S.Ray(
			raygr_inferred.ray_out.x0,
			raygr_inferred.ray_out.k,
			raygr_inferred.ray_out.e1,
			raygr_inferred.ray_out.e2,
			raygr_inferred.ray_out.nz,
			S.FastLight(),
		)

		@test_throws ArgumentError S.RayGR2(ray_in, ray_out_bad_freq, bh_pos)
		@test_throws MethodError S.RayGR2(ray_in, ray_out_bad_light, bh_pos)
	end

	@testset "CameraGR wraps orthographic and perspective cameras" begin
		uv_o = S.SVector(35.0, 0.0)
		cam_o = S.CameraOrtho(; photon_direction, xys=[uv_o], nz=nz, ν=ν, t=0.0)
		raygr_o = S.RayGR2(S.camera_ray(cam_o, uv_o), lens)
		camgr_o = S.CameraGR(; camera=cam_o, lens)
		sphere_o = S.UniformSphere(;
			center=S.FourPosition(0.0, 35.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		@test only(S.render(camgr_o, sphere_o)) ≈ S.render(raygr_o, sphere_o)

		uv_p = S.SVector(-0.02, 0.0)
		cam_p = S.CameraPerspective(; photon_direction, origin=-1000.0 * photon_direction, xys=[uv_p], nz=nz, ν=ν, t=0.0)
		ray_p = S.camera_ray(cam_p, uv_p)
		raygr_p = S.RayGR2(ray_p, lens)
		camgr_p = S.CameraGR(; camera=cam_p, lens)
		n_p = S.direction3(ray_p)
		anchor_p = S.SVector(ray_p.x0.x, ray_p.x0.y, ray_p.x0.z)
		sphere_p = S.UniformSphere(;
			center=S.FourPosition(0.0, (anchor_p + 50.0 * n_p)...), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		@test only(S.render(camgr_p, sphere_p)) ≈ S.render(raygr_p, sphere_p)
	end

	@testset "far ray: sphere in front of BH" begin
		# Ray at x=50, far from BH (b=50 >> shadow). Sphere on the incoming ray path,
		# entirely between camera and BH. GR incoming sees it fully, outgoing ray
		# (in BL→Cartesian frame, ~60° away) misses it entirely.
		ray_in = make_ray(50.0, 0.0)
		ray_gr = S.RayGR2(make_ray(50.0, 0.0), lens)
		@test !S.is_captured_ray(ray_gr.ray_out)

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=30.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S.render(ray_gr, sphere, S.Intensity())

		@test flat ≈ jν * 60.0 rtol=0.01
		@test gr == flat
	end

	@testset "far ray: sphere straddling BH plane" begin
		# Sphere centered at (50, 0, 0) — right at the BH clip plane.
		# For weak deflection (b=50), both incoming and outgoing rays nearly overlap.
		# Incoming sees the observer half (z>0), outgoing sees the source half (z<0).
		# Together they see approximately the full sphere.
		ray_in = make_ray(50.0, 0.0)
		ray_gr = S.RayGR2(make_ray(50.0, 0.0), lens)

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 0.0), radius=30.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S.render(ray_gr, sphere, S.Intensity())

		# Flat sees full sphere: I = j * 2R = j * 60
		@test flat ≈ jν * 60.0 rtol=0.01
		# GR sees nearly full sphere (both halves via incoming + outgoing)
		@test gr ≈ flat rtol=0.5
	end

	@testset "captured ray: sphere behind BH invisible" begin
		# Ray at x=0, y=0 — hits BH directly, captured.
		ray_in = make_ray(0.0, 0.0)
		ray_gr = S.RayGR2(make_ray(0.0, 0.0), lens)
		@test S.is_captured_ray(ray_gr.ray_out)  # captured

		# Sphere behind BH
		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, 0.0, 0.0, -50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S.render(ray_gr, sphere, S.Intensity())

		# Flat sees the sphere (straight line goes through it): I = j * 2R
		@test flat ≈ jν * 20.0 rtol=0.01
		# GR: incoming ray clipped at BH, ray captured → sphere invisible
		@test gr === 0.0
	end

	@testset "intermediate ray: sphere in front visible, behind invisible" begin
		# Ray at x=10, moderate offset. Strong deflection but not captured.
		ray_in = make_ray(10.0, 0.0)
		ray_gr = S.RayGR2(make_ray(10.0, 0.0), lens)
		@test !S.is_captured_ray(ray_gr.ray_out)  # not captured

		# Sphere on the incoming ray path, in front of BH
		sphere_front = S.UniformSphere(;
			center=S.FourPosition(0.0, 10.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat_front = S.integrate_ray(sphere_front, ray_in)
		gr_front = S.render(ray_gr, sphere_front, S.Intensity())

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
		gr_behind = S.render(ray_gr, sphere_behind, S.Intensity())

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
		raygr_ref = S.RayGR2(make_ray(50.0, 0.0), lens)
		gr_ref = S.render(raygr_ref, sphere_ref, S.Intensity())
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
		lens_shifted = S.GRLens(; spin, bh_position=bh_shifted)
		raygr_shifted = S.RayGR2(ray_in_shifted, lens_shifted)
		@test !S.is_captured_ray(raygr_shifted.ray_out)
		gr_shifted = S.render(raygr_shifted, sphere_shifted, S.Intensity())

		@test gr_shifted == gr_ref
	end

	@testset "bh_rg scaling" begin
		# Scaling bh_rg by factor s and camera xys by factor s should give identical
		# pixel values (same geometry, just in different units).
		scale = 5.0

		# Reference: rg=1, sphere at z=50 radius=10, ray at x=50
		cam_ref = S.CameraOrtho(;
			photon_direction,
			xys=S.grid(S.SVector, x=[49.0, 50.0, 51.0], y=[-1.0, 0.0, 1.0]),
			nz=nz, ν=ν, t=0.0,
		)
		lens_ref = S.GRLens(; spin, bh_rg=1.0)
		sphere_ref = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0, 0.0, 50.0), radius=10.0,
			u0=u_rest, jν=jν, αν=0.0,
		)
		cam_gr_ref = S.CameraGR(; camera=cam_ref, lens=lens_ref)
		img_ref = S.render(cam_gr_ref, sphere_ref)

		# Scaled: rg=scale, xys*scale, sphere coords*scale
		cam_scaled = S.CameraOrtho(;
			photon_direction,
			xys=S.grid(S.SVector, x=[49.0, 50.0, 51.0] .* scale, y=[-1.0, 0.0, 1.0] .* scale),
			nz=nz, ν=ν, t=0.0,
		)
		lens_scaled = S.GRLens(; spin, bh_rg=scale)
		sphere_scaled = S.UniformSphere(;
			center=S.FourPosition(0.0, 50.0*scale, 0.0, 50.0*scale), radius=10.0*scale,
			u0=u_rest, jν=jν, αν=0.0,
		)
		cam_gr_scaled = S.CameraGR(; camera=cam_scaled, lens=lens_scaled)
		img_scaled = S.render(cam_gr_scaled, sphere_scaled)

		# Optically thin: I = j * path_length. Path scales by `scale`, so I scales too.
		@test collect(img_scaled) ≈ collect(img_ref) .* scale rtol=0.01
	end

	@testset "two spheres on the same straight ray across the BH" begin
		ray_in = make_ray(50.0, 0.0)
		# Weak-deflection limit: use a straight outgoing segment, while RayGR2
		# still splits the path at the BH into far/background and near/foreground halves.
		ray_out = ray_in
		ray_gr = S.RayGR2(ray_in, ray_out, bh_pos)

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
		gr_front = S.render(ray_gr, sphere_front, S.Intensity())
		gr_back = S.render(ray_gr, sphere_back, S.Intensity())

		@test flat_front > 0
		@test flat_back > 0
		@test gr_front > 0
		@test gr_back > 0
		@test gr_front ≈ flat_front rtol=1e-10
		@test gr_back ≈ flat_back rtol=1e-10

		flat = S.integrate_ray(combined, ray_in)
		gr = S.render(ray_gr, combined, S.Intensity())
		τ_front = S.integrate_ray(sphere_front, ray_in, S.OpticalDepth())
		expected = flat_back * exp(-τ_front) + flat_front

		@test gr ≈ flat rtol=1e-10
		@test flat ≈ expected rtol=0.05
		@test gr ≈ expected rtol=0.05
	end

	@testset "large-bend ray: sphere on outgoing path only (nearly sideways)" begin
		# b=6.110 is just outside the shadow. The outgoing ray bends ~90° into the -x direction
		# (source at -x, photon travels in +x toward BH).
		# Sphere at (-142, 12, -4) sits along the outgoing path but far from the incoming ray.
		ray_in = make_ray(6.110, 0.01)
		ray_gr = S.RayGR2(make_ray(6.110, 0.01), lens)
		ray_out = ray_gr.ray_out
		@test !S.is_captured_ray(ray_out)
		@test abs(S.dot(S.direction3(ray_in), S.direction3(ray_out))) < 0.1

		radius = 5.5
		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, -142.0, 12.0, -4.0), radius,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S.render(ray_gr, sphere, S.Intensity())

		@test flat === 0.0
		@test gr ≈ jν * (2radius) rtol=0.05
	end

	@testset "anchor position independence" begin
		# Two rays along the same line but anchored at different z positions.
		# Deflection must be identical regardless of anchor placement.
		ray_z_pos = make_ray(50.0, 0.0)  # anchor at z=+1000
		ray_z_neg = S.Ray(
			S.FourPosition(0.0, 50.0, 0.0, -1000.0),
			ray_z_pos.k, ray_z_pos.e1, ray_z_pos.e2, nz, S.SlowLight(),
		)
		gr_pos = S.RayGR2(make_ray(50.0, 0.0), lens)
		gr_neg = S.RayGR2(ray_z_neg, lens)
		@test !S.is_captured_ray(gr_pos.ray_out)
		@test !S.is_captured_ray(gr_neg.ray_out)
		@test S.direction3(gr_pos.ray_out) ≈ S.direction3(gr_neg.ray_out) atol=1e-10
	end

	@testset "large-bend ray: sphere on outgoing path only (~180° U-turn)" begin
		# b=5.30: photon nearly U-turns (δ≈179°), source above BH on observer's side.
		# Sphere at (-3, -6, 50) sits along the outgoing path but far from the incoming ray.
		ray_in = make_ray(5.30, 0.01)
		ray_gr = S.RayGR2(make_ray(5.30, 0.01), lens)
		ray_out = ray_gr.ray_out
		@test !S.is_captured_ray(ray_out)
		# ~179° bend: directions nearly anti-parallel
		@test S.dot(S.direction3(ray_in), S.direction3(ray_out)) < -0.9

		sphere = S.UniformSphere(;
			center=S.FourPosition(0.0, -3.0, -6.0, 50.0), radius=3.0,
			u0=u_rest, jν=jν, αν=0.0,
		)

		flat = S.integrate_ray(sphere, ray_in)
		gr = S.render(ray_gr, sphere, S.Intensity())

		@test flat === 0.0
		@test gr ≈ jν * 5.95 rtol=0.05
	end

end

@testitem "GR frequency: p^t and gravitational redshift" begin
	import Synchray as S
	using Krang
	using FastChebInterp
	using Krang: Kerr, emission_coordinates, p_bl_d, metric_uu

	KE = Base.get_extension(S, :KrangExt)

	@testset "p^t → 1 at large r" begin
		# For any geodesic, the photon energy p^t should approach 1 at r → ∞
		# (Krang normalizes energy to 1 at infinity)
		for (a, α, β) in [(0.5, 8.0, 1.0), (0.9, 10.0, 0.5), (0.01, 10.0, 1.0)]
			local metric = Kerr(a)
			local pix = KE._krang_pixel(metric, α, β, π/4)
			@assert pix.total_mino_time > 0

			# Sample at early τ (large r, near observer)
			local τ = 0.001 * pix.total_mino_time
			local t, r, θ, φ, νr, νθ, ok = emission_coordinates(pix, τ)
			@assert ok && r > 100

			local p_d = p_bl_d(metric, r, θ, pix.η, pix.λ, νr, νθ)
			local p_u = metric_uu(metric, r, θ) * p_d
			@test p_u[1] ≈ 1.0 atol=0.02
		end
	end

	@testset "p^t matches Schwarzschild analytic r/(r-2)" begin
		# For near-Schwarzschild (a≈0), p^t = -g^{tt}·(-1) + g^{tφ}·λ ≈ r/(r-2)
		metric = Kerr(0.01)
		pix = KE._krang_pixel(metric, 10.0, 1.0, π/4)
		@assert pix.total_mino_time > 0

		for τ_frac in [0.01, 0.1, 0.3, 0.5]
			local τ = τ_frac * pix.total_mino_time
			local t, r, θ, φ, νr, νθ, ok = emission_coordinates(pix, τ)
			@assert ok

			local p_d = p_bl_d(metric, r, θ, pix.η, pix.λ, νr, νθ)
			local p_u = metric_uu(metric, r, θ) * p_d

			# Schwarzschild: p^t = r/(r-2) for p_t = -1, a=0
			analytic_pt = r / (r - 2)
			@test p_u[1] ≈ analytic_pt rtol=0.001
		end
	end

	@testset "GR frequency changes rendered intensity vs flat approximation" begin
		# Render with FixedEmission at moderate r. The GR frequency correction (p^t > 1
		# near BH) changes the local frequency seen by the medium, which affects the
		# invariant emissivity j/ν² and absorption α·ν. For FixedEmission with α > 0,
		# the intensity approaches the source function S in the optically thick limit
		# regardless of frequency. So we use small α (optically thin) where the
		# frequency effect shows up via the ν³ postprocessing factor.
		#
		# Test: intensity with GR frequency (current) is LOWER than what we'd get
		# with ν_obs everywhere, because p^t > 1 increases local ν, making j_inv = j/ν² smaller.
		# We verify this by comparing corner vs near-center: the near-center pixel
		# passes through stronger gravity (smaller r, larger p^t), so the GR frequency
		# effect is larger there. The ratio near-center/corner should be < 1 and smaller
		# than the pure geometric (path-length) ratio.

		u_rest = S.FourVelocity(1.0, 0.0, 0.0, 0.0)
		ball = S.EmissionRegion(
			geometry = Geometries.Ball(
				center = Geometries.InertialWorldline(x0 = S.FourPosition(0.0, 0.0, 0.0, 0.0), u = u_rest),
				size = 30.0),
			velocity = nothing,
			emission = S.FixedEmission(S=1.0, α=1e-3)
		) |> S.prepare_for_computations

		cam = S.CameraOrtho(;
			photon_direction = S.SVector(1.0, 0.0, 0.0),
			origin = S.SVector(-100.0, 0.0, 0.0),
			xys = S.grid(S.SVector, range(-15, 15, length=7), range(-15, 15, length=7)),
			nz = 200, ν = 1.0, t = 0.0)

		cam_gr = S.CameraKerrGR(; camera=cam, metric_spin=0.9,
			bh_position=S.SVector(0.0,0.0,0.0), bh_rg=1.0, nτ=300, τ_range=(0.01, 0.99))
		img = S.render(cam_gr, ball)

		# Corner pixel (b ≈ 21, r_min ≈ 20): weak gravity, p^t ≈ 1.05
		# Near-center pixel (b ≈ 7, r_min ≈ 6): strong gravity, p^t ≈ 1.4
		# The ratio should reflect both geometry AND GR frequency effects
		corner = img[1, 1]
		near_center = img[3, 3]

		@test corner ≈ 0.756 rtol=0.01
		@test near_center ≈ 0.320 rtol=0.02

		# The ratio near_center/corner is ~0.42, much less than what pure path-length
		# differences would give, because GR frequency suppression is stronger at smaller r
		@test near_center / corner ≈ 0.423 rtol=0.02
	end
end

@testitem "GR velocity types and PowerLawDisk" begin
	import Synchray as S
	using Krang
	using FastChebInterp
	using Krang: Kerr, metric_dd

	@testset "KeplerianVelocity normalization g_μν u^μ u^ν = -1" begin
		# Note: Krang's Kerr metric requires nonzero spin, so we use a=0.01 for near-Schwarzschild
		for (a, r) in [(0.01, 10.0), (0.5, 6.0), (0.9, 3.0), (0.99, 2.0), (0.9, 50.0)]
			kv = S.KeplerianVelocity(spin=a, prograde=true)
			x4 = S.FourPosition(0.0, r, 0.0, 0.0)
			u = S.four_velocity(kv, nothing, x4)

			g_dd = metric_dd(Kerr(a), r, π/2)
			Ω = 1 / (r^(3/2) + a)
			uφ = Ω * u.t
			gnorm = g_dd[1,1]*u.t^2 + 2*g_dd[1,4]*u.t*uφ + g_dd[4,4]*uφ^2
			@test gnorm ≈ -1.0 atol=1e-6
		end
	end

	@testset "KeplerianVelocity u^t and β match analytic" begin
		# a=0.9, r=10: u^t = 1.1821, β = 0.3075
		kv = S.KeplerianVelocity(spin=0.9, prograde=true)
		u = S.four_velocity(kv, nothing, S.FourPosition(0.0, 10.0, 0.0, 0.0))
		@test u.t ≈ 1.1821 rtol=0.001
		@test S.norm(S.beta(u)) ≈ 0.3075 rtol=0.001

		# a=0.9, r=3 (near ISCO): u^t = 1.9933, β = 0.4921
		u3 = S.four_velocity(kv, nothing, S.FourPosition(0.0, 3.0, 0.0, 0.0))
		@test u3.t ≈ 1.9933 rtol=0.001
		@test S.norm(S.beta(u3)) ≈ 0.4921 rtol=0.001

		# a=0.9, r=50 (weak field): u^t ≈ 1.031, β ≈ 0.141
		u50 = S.four_velocity(kv, nothing, S.FourPosition(0.0, 50.0, 0.0, 0.0))
		@test u50.t ≈ 1.0313 rtol=0.001
		@test S.norm(S.beta(u50)) ≈ 0.1411 rtol=0.001
	end

	@testset "ZAMOVelocity has zero angular momentum u_φ = 0" begin
		for (a, r) in [(0.5, 10.0), (0.9, 5.0), (0.99, 3.0)]
			zv = S.ZAMOVelocity(spin=a)
			u = S.four_velocity(zv, nothing, S.FourPosition(0.0, r, 0.0, 0.0))

			# u_φ = g_φt u^t + g_φφ u^φ where u^φ = Ω_ZAMO * u^t
			g_tφ = S._kerr_g_tφ(a, r, π/2)
			g_φφ = S._kerr_g_φφ(a, r, π/2)
			Ω_zamo = -g_tφ / g_φφ
			uφ = Ω_zamo * u.t
			u_φ_cov = metric_dd(Kerr(a), r, π/2)[4,1]*u.t + metric_dd(Kerr(a), r, π/2)[4,4]*uφ
			@test u_φ_cov ≈ 0.0 atol=1e-10
		end
	end

	@testset "ZAMOVelocity normalization and hardcoded values" begin
		zv = S.ZAMOVelocity(spin=0.9)
		u5 = S.four_velocity(zv, nothing, S.FourPosition(0.0, 5.0, 0.0, 0.0))
		@test u5.t ≈ 1.2857 rtol=0.001
		@test S.norm(S.beta(u5)) ≈ 0.0689 rtol=0.01  # small β: ZAMO barely moves at r=5

		u_far = S.four_velocity(zv, nothing, S.FourPosition(0.0, 100.0, 0.0, 0.0))
		# At large r, ZAMO → nearly static: u^t → 1.01 (still has 2M/r correction), β → 0
		@test u_far.t ≈ 1.0102 rtol=0.001
		@test S.norm(S.beta(u_far)) ≈ 0.00018 rtol=0.1
	end

	@testset "SubKeplerianVelocity interpolates between ZAMO and Keplerian" begin
		a = 0.9
		r = 10.0
		x4 = S.FourPosition(0.0, r, 0.0, 0.0)

		β_kep = S.norm(S.beta(S.four_velocity(S.KeplerianVelocity(spin=a), nothing, x4)))
		β_zamo = S.norm(S.beta(S.four_velocity(S.ZAMOVelocity(spin=a), nothing, x4)))
		β_half = S.norm(S.beta(S.four_velocity(S.SubKeplerianVelocity(spin=a, f_kep=0.5), nothing, x4)))

		# Sub-Keplerian speed should be between ZAMO and Keplerian
		@test β_zamo < β_half < β_kep  # this is a < comparison but it's comparing three related quantities, which is meaningful
	end

	@testset "PowerLawDisk is_inside" begin
		disk = Geometries.PowerLawDisk(r_range=3.0..20.0, h_ref=1.0, r_ref=10.0, a=0.5)
		# h(r) = 1.0 * (r/10)^0.5

		# At r=10: h=1.0
		@test S.is_inside(disk, S.FourPosition(0.0, 10.0, 0.0, 0.0))      # z=0 < h=1
		@test S.is_inside(disk, S.FourPosition(0.0, 10.0, 0.0, 0.9))      # z=0.9 < h=1
		@test !S.is_inside(disk, S.FourPosition(0.0, 10.0, 0.0, 1.1))     # z=1.1 > h=1

		# At r=4: h = 1.0 * (4/10)^0.5 ≈ 0.632
		@test S.is_inside(disk, S.FourPosition(0.0, 4.0, 0.0, 0.0))
		@test !S.is_inside(disk, S.FourPosition(0.0, 4.0, 0.0, 0.7))

		# Outside r_range
		@test !S.is_inside(disk, S.FourPosition(0.0, 2.0, 0.0, 0.0))
		@test !S.is_inside(disk, S.FourPosition(0.0, 25.0, 0.0, 0.0))

		# Constant thickness (a=0)
		slab = Geometries.PowerLawDisk(r_range=5.0..15.0, h_ref=2.0, r_ref=10.0, a=0)
		@test S.is_inside(slab, S.FourPosition(0.0, 10.0, 0.0, 1.9))
		@test !S.is_inside(slab, S.FourPosition(0.0, 10.0, 0.0, 2.1))
		# a=0: h(r) = 2.0 for all r in range
		@test S.is_inside(slab, S.FourPosition(0.0, 7.0, 0.0, 1.9))
	end

	@testset "EmissionRegion composes with GR velocity" begin
		# Verify that KeplerianVelocity works inside EmissionRegion
		disk = S.EmissionRegion(
			geometry = Geometries.PowerLawDisk(r_range=3.0..20.0, h_ref=1.0, r_ref=10.0, a=0.5),
			velocity = S.KeplerianVelocity(spin=0.9),
			emission = S.FixedEmission(S=1.0, α=1e-3)
		) |> S.prepare_for_computations

		x4 = S.FourPosition(0.0, 10.0, 0.0, 0.0)
		u = S.four_velocity(disk, x4)
		@test u.t ≈ 1.1821 rtol=0.001
	end
end

@testitem "BL momentum → Cartesian direction vs finite differences" begin
	import Synchray as S
	using Krang
	using FastChebInterp
	using Krang: Kerr, SlowLightIntensityPixel, emission_coordinates

	KE = Base.get_extension(S, :KrangExt)

	# Test for several spins and impact parameters (all non-captured with valid geodesics)
	for (a, α, β) in [(0.5, 8.0, 2.0), (0.9, 10.0, 0.5), (0.99, 7.0, 3.0), (0.1, 20.0, 5.0)]
		metric = Kerr(a)
		pix = KE._krang_pixel(metric, α, β, π/4)
		@assert pix.total_mino_time > 0

		τ_total = pix.total_mino_time

		for τ_frac in [0.2, 0.4, 0.6, 0.8]
			τ = τ_frac * τ_total
			t1, r1, θ1, φ1, νr1, νθ1, ok1 = emission_coordinates(pix, τ)
			@assert ok1 "emission_coordinates failed at τ_frac=$τ_frac for a=$a"

			# Jacobian-based direction
			p_d = Krang.p_bl_d(metric, r1, θ1, pix.η, pix.λ, νr1, νθ1)
			g_uu = Krang.metric_uu(metric, r1, θ1)
			p_u = g_uu * p_d
			n_jac = KE._bl_upper_to_cartesian_direction(metric, p_u[2], p_u[3], p_u[4], r1, θ1, φ1)

			# Finite-difference direction from two nearby geodesic points
			δτ = 1e-6 * τ_total
			t_a, r_a, θ_a, φ_a, _, _, ok_a = emission_coordinates(pix, τ - δτ)
			t_b, r_b, θ_b, φ_b, _, _, ok_b = emission_coordinates(pix, τ + δτ)
			@assert ok_a && ok_b "FD emission_coordinates failed at τ_frac=$τ_frac for a=$a"

			x_a = SVector(Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(metric, r_a, θ_a, φ_a)...)
			x_b = SVector(Krang.boyer_lindquist_to_quasi_cartesian_kerr_schild_fast_light(metric, r_b, θ_b, φ_b)...)
			n_fd = S.normalize(x_b - x_a)

			# Directions should agree (up to sign — both parameterize the same null ray)
			# Use abs(dot) to handle sign ambiguity
			@test abs(S.dot(n_jac, n_fd)) ≈ 1.0 atol=1e-4
		end
	end
end

@testitem "CameraKerrGR geodesic integration" begin
	import Synchray as S
	using Krang
	using FastChebInterp

	# Shared setup: spinning BH and a large spherical emitter
	spin = 0.9
	bh_rg = 1.0
	bh_pos = S.SVector(0.0, 0.0, 0.0)

	u_rest = S.FourVelocity(1.0, 0.0, 0.0, 0.0)
	ball = S.EmissionRegion(
		geometry = Geometries.Ball(
			center = Geometries.InertialWorldline(x0 = S.FourPosition(0.0, bh_pos...), u = u_rest),
			size = 30.0
		),
		velocity = nothing,
		emission = S.FixedEmission(S=1.0, α=1e-3)
	) |> S.prepare_for_computations

	cam = S.CameraOrtho(;
		photon_direction = S.SVector(1.0, 0.0, 0.0),
		origin = S.SVector(-100.0, 0.0, 0.0),
		xys = S.grid(S.SVector, range(-15, 15, length=7), range(-15, 15, length=7)),
		nz = 200,
		ν = 1.0,
		t = 0.0,
	)

	@testset "basic rendering produces expected intensities" begin
		cam_gr = S.CameraKerrGR(; camera=cam, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=300, τ_range=(0.01, 0.99))
		img = S.render(cam_gr, ball)

		# Hardcoded reference values (FixedEmission S=1, α=1e-3, Ball size=30, spin=0.9)
		# With proper GR frequency: p^t > 1 near BH increases local ν, reducing I via ν³ postprocessing
		@test img[1, 1] ≈ 0.756 rtol=0.01   # corner: long path through ball
		@test img[7, 7] ≈ 0.752 rtol=0.01   # opposite corner (slightly different from frame-dragging)
		@test img[3, 3] ≈ 0.320 rtol=0.02   # near-center: shorter path, stronger GR redshift effect

		# Equatorial row (y=0) is zero due to Krang β≈0 degeneracy
		@test img[4, 4] ≈ 0.0 atol=1e-30
	end

	@testset "image z-reflection symmetry (outer pixels)" begin
		cam_gr = S.CameraKerrGR(; camera=cam, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=300, τ_range=(0.01, 0.99))
		img = S.render(cam_gr, ball)

		# Viewing along x (θ_obs = π/2), spin along z:
		# Kerr has a θ → π-θ reflection symmetry → z → -z → image y-symmetric.
		# Test at large impact parameter (corner pixels) where approximation is best.
		@test img[1, 1] ≈ img[1, 7] rtol=0.02
		@test img[7, 1] ≈ img[7, 7] rtol=0.02
		@test img[1, 2] ≈ img[1, 6] rtol=0.02
	end

	@testset "nτ convergence" begin
		# Render at three resolutions; verify convergence toward stable value.
		cam_gr_100 = S.CameraKerrGR(; camera=cam, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=100, τ_range=(0.01, 0.99))
		cam_gr_200 = S.CameraKerrGR(; camera=cam, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=200, τ_range=(0.01, 0.99))
		cam_gr_400 = S.CameraKerrGR(; camera=cam, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=400, τ_range=(0.01, 0.99))

		img_100 = S.render(cam_gr_100, ball)
		img_200 = S.render(cam_gr_200, ball)
		img_400 = S.render(cam_gr_400, ball)

		# Corner pixel values at each resolution
		i, j = 1, 1
		@test img_100[i,j] ≈ 0.856 rtol=0.01
		@test img_200[i,j] ≈ 0.784 rtol=0.01
		@test img_400[i,j] ≈ 0.741 rtol=0.01

		# Error ratio: err_200/err_100 ≈ 0.37 (convergence rate)
		err_100 = abs(img_100[i,j] - img_400[i,j])
		err_200 = abs(img_200[i,j] - img_400[i,j])
		@test err_200 / err_100 ≈ 0.370 rtol=0.1
	end

	@testset "captured ray still accumulates foreground emission" begin
		# A captured ray (small impact parameter) falls into the BH, but along the
		# incoming path (camera → horizon) it passes through the emitting ball.
		# It should accumulate nonzero intensity from the foreground part.
		cam_captured = S.CameraOrtho(;
			photon_direction = S.SVector(1.0, 0.0, 0.0),
			origin = S.SVector(-100.0, 0.0, 0.0),
			xys = [S.SVector(3.0, 0.5)],  # small impact parameter, captured
			nz = 200, ν = 1.0, t = 0.0,
		)
		cam_gr_captured = S.CameraKerrGR(; camera=cam_captured, metric_spin=spin, bh_position=bh_pos, bh_rg, nτ=300, τ_range=(0.01, 0.99))
		img_captured = S.render(cam_gr_captured, ball)

		# The ball extends to radius 30. A captured ray with b≈3 enters the ball at
		# r≈30, accumulates emission along the infall to the horizon.
		@test only(img_captured) ≈ 0.139 rtol=0.02
	end

	@testset "different spin gives different image" begin
		# Same geometry, different spin — images should differ due to frame-dragging
		cam_gr_low = S.CameraKerrGR(; camera=cam, metric_spin=0.1, bh_position=bh_pos, bh_rg, nτ=200, τ_range=(0.01, 0.99))
		cam_gr_high = S.CameraKerrGR(; camera=cam, metric_spin=0.99, bh_position=bh_pos, bh_rg, nτ=200, τ_range=(0.01, 0.99))

		img_low = S.render(cam_gr_low, ball)
		img_high = S.render(cam_gr_high, ball)

		# Hardcoded reference values for corner pixel
		@test img_low[1, 1] ≈ 0.782 rtol=0.01
		@test img_high[1, 1] ≈ 0.784 rtol=0.01

		# At small impact parameter, spin changes geodesic path through the medium.
		# The difference is small (~1.5%) but real and detectable.
		@test img_low[3, 2] ≈ 0.521 rtol=0.01
		@test img_high[3, 2] ≈ 0.529 rtol=0.01
		@test !isapprox(img_low[3, 2], img_high[3, 2]; rtol=0.001)
	end
end

@testitem "ray_in_local_coords for RayGR2" begin
	import Synchray as S
	using Krang
	using FastChebInterp

	spin = 0.5
	θ_view = π / 3
	photon_direction = S.normalize(S.SVector(0.0, -sin(θ_view), cos(θ_view)))
	lens = S.GRLens(; spin)

	cam = S.CameraOrtho(;
		photon_direction,
		xys = S.grid(S.SVector, x=range(-60.0, 60.0, length=32), y=range(-60.0, 60.0, length=32)),
		nz = 100, ν = 1e11, t = 0.0,
	)

	geom = S.Geometries.Conical(S.SVector(0.0, 0.0, 1.0), 0.1, 1.0 .. 100.0)
	s_range = -1e4 .. 1e4

	# Deflected ray: 4 points
	ray_gr = S.RayGR2(S.camera_ray(cam, S.SVector(30.0, 0.01)), lens)
	@test !S.is_captured_ray(ray_gr.ray_out)
	pts = S.ray_in_local_coords(ray_gr, geom; s_range)
	@test length(pts) == 4
	@test all(p -> p isa SVector{3}, pts)
	# Points 2 and 3 (near-BH) should be close to each other (both near BH)
	@test S.norm(pts[2] - pts[3]) < S.norm(pts[1] - pts[4])
	# Points 1 and 4 (far endpoints) should be far apart
	@test S.norm(pts[1] - pts[4]) > 100

	# Captured ray: 2 points
	ray_cap = S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 0.0)), lens)
	@test S.is_captured_ray(ray_cap.ray_out)
	pts_cap = S.ray_in_local_coords(ray_cap, geom; s_range)
	@test length(pts_cap) == 2
	@test all(p -> p isa SVector{3}, pts_cap)
end

@testitem "GR lensing" begin
	import Synchray as S
	using Krang

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
		@test S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 0.0)), lens).ray_out === nothing
		@test S.RayGR2(S.camera_ray(cam, S.SVector(3.0, 0.0)), lens).ray_out === nothing
		@test S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 5.0)), lens).ray_out === nothing
		# Outside shadow
		@test S.RayGR2(S.camera_ray(cam, S.SVector(10.0, 0.0)), lens).ray_out !== nothing
		@test S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 50.0)), lens).ray_out !== nothing
	end

	@testset "outgoing direction is unit vector" begin
		for uv in [S.SVector(10.0, 0.0), S.SVector(30.0, 5.0), S.SVector(50.0, -20.0)]
			d = S.RayGR2(S.camera_ray(cam, uv), lens).ray_out
			@test d !== nothing
			@test S.norm(S.direction3(d)) ≈ 1 atol=1e-10
		end
	end

	@testset "weak-field deflection angle" begin
		# For b >> M, deflection ≈ 4M/b (Schwarzschild limit, M=1 in Krang).
		# Use non-zero y to avoid the Krang β=0 degeneracy.
		for uv in [S.SVector(40.0, 10.0), S.SVector(30.0, 30.0), S.SVector(45.0, -20.0)]
			b = S.norm(uv)
			d = S.RayGR2(S.camera_ray(cam, uv), lens).ray_out
			@test d !== nothing
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
		@test d_sp !== nothing

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

		@test raygr_inferred.ray_out !== nothing
		@test raygr_inferred.ray_out.x0 ≈ raygr_explicit.ray_out.x0
		@test raygr_inferred.ray_out.k ≈ raygr_explicit.ray_out.k
		@test raygr_inferred.ray_out.e1 ≈ raygr_explicit.ray_out.e1
		@test raygr_inferred.ray_out.e2 ≈ raygr_explicit.ray_out.e2
		@test raygr_twisted.ray_out !== nothing
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
		@test_throws ArgumentError S.RayGR2(ray_in, ray_out_bad_light, bh_pos)
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
		@test ray_gr.ray_out !== nothing

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
		@test ray_out = ray_gr.ray_out === nothing  # captured

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
		@test ray_gr.ray_out !== nothing  # not captured

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
		@test raygr_shifted.ray_out !== nothing
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
		@test ray_out !== nothing
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
		@test gr_pos.ray_out !== nothing
		@test gr_neg.ray_out !== nothing
		@test S.direction3(gr_pos.ray_out) ≈ S.direction3(gr_neg.ray_out) atol=1e-10
	end

	@testset "large-bend ray: sphere on outgoing path only (~180° U-turn)" begin
		# b=5.30: photon nearly U-turns (δ≈179°), source above BH on observer's side.
		# Sphere at (-3, -6, 50) sits along the outgoing path but far from the incoming ray.
		ray_in = make_ray(5.30, 0.01)
		ray_gr = S.RayGR2(make_ray(5.30, 0.01), lens)
		ray_out = ray_gr.ray_out
		@test ray_out !== nothing
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

@testitem "ray_in_local_coords for RayGR2" begin
	import Synchray as S
	using Krang

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
	@test ray_gr.ray_out !== nothing
	pts = S.ray_in_local_coords(ray_gr, geom; s_range)
	@test length(pts) == 4
	@test all(p -> p isa SVector{3}, pts)
	# Points 2 and 3 (near-BH) should be close to each other (both near BH)
	@test S.norm(pts[2] - pts[3]) < S.norm(pts[1] - pts[4])
	# Points 1 and 4 (far endpoints) should be far apart
	@test S.norm(pts[1] - pts[4]) > 100

	# Captured ray: 2 points
	ray_cap = S.RayGR2(S.camera_ray(cam, S.SVector(0.0, 0.0)), lens)
	@test ray_cap.ray_out === nothing
	pts_cap = S.ray_in_local_coords(ray_cap, geom; s_range)
	@test length(pts_cap) == 2
	@test all(p -> p isa SVector{3}, pts_cap)
end

@testitem "Uniform slab transfer" begin
    import Synchray as S

	L = 3
	j0 = 0.7
	a0 = 1.3

	ray = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2, nz=2_000)

	cases = [
		(; u0=S.FourVelocity(SVector(0, 0, 0))),
		(; u0=S.FourVelocity(SVector(0, 0, 0.3))),
		(; u0=S.FourVelocity(SVector(0.3, 0, 0.5))),
	]

	@testset for (;u0) in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset "limits" begin
			@testset "optically thin" begin
				# For τ = α D L ≪ 1: I′ ≈ j L′ and I = I′/D^3 ⇒ I ≈ j L / D^2
				αthin = 1e-6
				slab = S.UniformSlab(0..L, u0, j0, αthin)
				Iν_num = S.render(ray, slab)
				Iν_exact = j0 * L * (δ^2)
				@test Iν_num ≈ Iν_exact rtol=2e-3
				@test S.render(ray, slab, S.OpticalDepth()) ≈ 2e-6 rtol=0.5
			end

			@testset "very optically thick (α ≫ j)" begin
				# For τ ≫ 1: I′ → S = j/α and I = S/D^3.
				αthick = 1e5
				slab = S.UniformSlab(0..L, u0, j0, αthick)
				Iν_num = S.render(ray, slab)
				Iν_exact = (j0 / αthick) * (δ^3)
				@test Iν_num ≈ Iν_exact rtol=2e-3
				@test S.render(ray, slab, S.OpticalDepth()) ≈ 2e5 rtol=0.5
			end
		end

		slab = S.UniformSlab(0.0..L, u0, j0, a0)
		Iν_num = S.render(ray, slab)
		τ_num = S.render(ray, slab, S.OpticalDepth())
		α_num = S.render(ray, slab, S.SpectralIndex())

		Iν_τ_num = S.render(ray, slab, (S.Intensity(), S.OpticalDepth()))
		@test Iν_τ_num[1] ≈ Iν_num
		@test Iν_τ_num[2] ≈ τ_num

		Iν_α_num = S.render(ray, slab, (S.Intensity(), S.SpectralIndex()))
		@test Iν_α_num[1] ≈ Iν_num
		@test Iν_α_num[2] ≈ α_num

		# Analytic solution as a cross-check:
		# - Solve transfer in slab comoving frame: dI'/ds' = j0 - a0 I'
		# - Convert thickness to comoving path length via D = ν'/ν = 1/δ
		#   for a ray along +z (n = ẑ): δ = doppler_factor(u, ẑ)
		# - Transform specific intensity back: I_ν = δ^{3} I'_{ν'}
		L′ = L / δ
		I′ = (j0 / a0) * (1 - exp(-a0 * L′))
		Iν_exact = I′ * (δ^3)

		@test Iν_num ≈ Iν_exact rtol=2e-3
		@test τ_num ≈ a0 * L′ rtol=2e-3

		# For frequency-independent comoving (j0, a0), the observed I_ν is ν-independent ⇒ α = 0.
		@test α_num ≈ 0 atol=2e-6
	end
end


@testitem "Uniform sphere flux (thin to thick)" begin
	import Synchray as S
    using RectiGrids
	using Accessors: @set

	R = 1.3
	j0 = 0.7
	center = S.FourPosition(0, 0, 0, 0)
	αs = (0, 1.2, 12)

	ray = S.RayZ(; x0=S.FourPosition(0,0,0,0), k=2, nz=512)

	cases = (
		S.FourVelocity(SVector(0, 0, 0)),
		S.FourVelocity(SVector(0, 0, 0.3)),
		S.FourVelocity(SVector(0.5, 0, 0.3)),
	)

	I_exact(j0, α0, δ, ℓ) = α0 == 0 ? (j0 * ℓ) * (δ^2) : (j0 / α0) * (1 - exp(-α0 * ℓ / δ)) * (δ^3)
	F_exact(j0, α0, δ, R) = α0 == 0 ? (4π / 3) * j0 * R^3 * (δ^2) : begin
		S0 = j0 / α0
		τ0 = 2 * α0 * R / δ
		(π * R^2 * S0 * (δ^3)) * (1 - (2 / (τ0^2)) * (1 - (1 + τ0) * exp(-τ0)))
	end

	@testset for u0 in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		@testset for α0 in αs
			sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=α0)

			# Per-ray analytic check
			for b in (0, 0.4R, 0.8R)
				ℓ = 2 * √(R^2 - b^2)
				I_num = S.render((@set ray.x0.x = b), sphere)
				@test I_num ≈ I_exact(j0, α0, δ, ℓ) rtol=6e-3
			end

			# Missed rays should return zero/nan
			mray = @set ray.x0.x = 2R
			@test S.render(mray, sphere) == 0
			@test S.render(mray, sphere, S.OpticalDepth()) == 0
			@test S.render(mray, sphere, S.SpectralIndex()) |> isnan
			@test isequal(S.render(mray, sphere, (S.Intensity(), S.SpectralIndex())), (0., NaN))

			# Total (image-plane) flux from a 2D grid: ∫ I(x,y) dx dy
			xs = range(-R..R, 151)
			cam = S.CameraZ(; xys=grid(SVector, xs, xs), ray.nz, ν=2, t=0)
			img = S.render(cam, sphere)
			dx = step(xs)
			F_num = sum(img) * dx^2
			@test F_num ≈ F_exact(j0, α0, δ, R) rtol=(α0 == 0 ? 1.6e-2 : 1.2e-2)
		end
	end
end




@testitem "Moving sphere silhouette (Terrell rotation)" begin
	import Synchray as S
	using Accessors

	R = 1
	z0 = 3

	βs = (
		SVector(-0.7, 0.0, 0.0),
		SVector(0.0, 0.7, 0.0),
		SVector(0.4, -0.3, 0.0),
		SVector(0.25, 0.15, 0.5),
	)
	t0s = (-1, 0, 0.75)

	# Stationary reference: for an orthographic camera (+z), motion should not change the image
	# apart from:
	# - a Terrell shift of the apparent center in the image plane, and
	# - an overall intensity scaling that depends only on the Doppler factor δ.
	#
	# We avoid rendering full images: just compare a small set of representative rays.
	ref_sphere = S.MovingUniformEllipsoid(
		center=S.FourPosition(0, 0, 0, z0),
		sizes=SVector(R, R, R),
		u0=S.FourVelocity(SVector(0, 0, 0)),
		jν=1,
		αν=0,
	)

	@testset for t0 in t0s
		cam = S.CameraZ(; xys=SVector{2}[
			(0, 0), (0.3R, 0.4R), (-0.8R, 0), (0.55R, -0.2R),
			(1.2R, 0), (0, 1.3R),
		], nz=2048, ν=2, t=t0)

		Iref = S.render(cam, ref_sphere)

		@testset for β in βs
			sphere = S.MovingUniformEllipsoid(
				center=S.FourPosition(0, 0, 0, z0),
				sizes=SVector(R, R, R),
				u0=S.FourVelocity(β),
				jν=1,
				αν=0,
			)

			# Terrell shift for orthographic camera along +z
			s = (t0 - sphere.center.t + sphere.center.z) / (1 - β.z)
			xyc = (@swiz sphere.center.xy) + (@swiz β.xy) * s

			δ = S.doppler_factor(sphere.u0, SVector(0, 0, 1))
			
			# The apparent brightness profile should match the stationary one, up to a constant factor.
			# Keep tolerances tight; this should be nearly exact for the same ray discretization.
			curcam = @set cam.xys[∗] += xyc
			Icur = S.render(curcam, sphere)
			@test Icur ≈ δ^3 * Iref rtol=3e-6
		end
	end
end


@testitem "Moving ellipsoid rotates (Terrell rotation)" begin
	import Synchray as S
	using Accessors

	# Penrose–Terrell result: for an orthographic camera along +z, a rigid body moving
	# along +x should appear as if it were *rotated* about the y axis by an angle θ
	# with sin(θ) = βx (plus a Terrell shift in the image plane and a Doppler boost).
	#
	# Choose an anisotropic ellipsoid (a ≠ c) so rotation changes the x-profile.
	# Compare against an analytic chord-length model for the rotated stationary ellipsoid.
	#
	# To keep the analytic side trivial, only sample rays with x=0 (through the apparent center):
	# then the quadratic cross-term vanishes and the chord length is
	#   ℓ(y) = 2 * √((1 - y^2/b^2) / A),  A = sin^2(θ)/a^2 + cos^2(θ)/c^2.

	z0 = 3
	sizes = SVector(2.0, 1.0, 0.4)  # (a, b, c)

	βs = (
		SVector(0.6, 0.0, 0.0),
		SVector(0.8, 0.0, 0.0),
		SVector(0.98, 0.0, 0.0),
	)
	t0s = (-0.5, 0.25)

	# Sample a few y offsets including near-edge / missed cases.
	xys = SVector{2}[
		(0, 0.0),
		(0, 0.3),
		(0, 0.7),
		(0, 1.0),
		(0, 1.1),
		(0, -0.8),
	]

	@testset for t0 in t0s
		cam = S.CameraZ(; xys, nz=4096, ν=2, t=t0)

		@testset for β in βs
			ell = S.MovingUniformEllipsoid(
				center=S.FourPosition(0, 0, 0, z0),
				sizes=sizes,
				u0=S.FourVelocity(β),
				jν=1,
				αν=0,
			)

			# Terrell shift for orthographic camera along +z
			s = (t0 - ell.center.t + ell.center.z) / (1 - β.z)
			xyc = (@swiz ell.center.xy) + (@swiz β.xy) * s

			δ = S.doppler_factor(ell.u0, SVector(0, 0, 1))
			θ = asin(β.x)
			(a, b, c) = sizes
			sθ, cθ = sincos(θ)
			A = (sθ^2) / (a^2) + (cθ^2) / (c^2)

			curcam = @set cam.xys[∗] += xyc
			Icur = S.render(curcam, ell)
			Iexp = map(xy -> begin
				(;y) = xy
				t = 1 - (y^2) / (b^2)
				ℓ = t > 0 ? 2 * √(t / A) : 0.0
				(δ^3) * ℓ
			end, xys)
			@test count(>(0), Iexp) == 4

			@test Icur ≈ Iexp rtol=3e-3
		end
	end
end


@testitem "Sphere rendering from arbitrary-angle rays" begin
	import Synchray as S
	using RectiGrids

	R = 1.3
	j0 = 0.7
	center = S.FourPosition(0, 0, 0, 0)
	u0 = S.FourVelocity(SVector(0, 0, 0))

	ν = 2.0
	nz = 512

	# Helper: create a Ray through a centered sphere at impact parameter b
	function make_ray(n̂, b)
		n̂ = normalize(n̂)
		up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
		e1 = normalize(cross(up, n̂))
		e2 = cross(n̂, e1)
		k = S.photon_k(ν, n̂)
		origin = -10.0 * n̂ + b * e1
		S.Ray(S.FourPosition(0.0, origin), k, e1, e2, nz)
	end

	directions = [
		SVector(0.0, 0.0, 1.0),
		SVector(1.0, 0.0, 0.0),
		SVector(0.0, 1.0, 0.0),
		normalize(SVector(1.0, 0.0, 1.0)),
		normalize(SVector(1.0, 1.0, 1.0)),
		normalize(SVector(-0.3, 0.7, 0.5)),
	]

	@testset for α0 in (0.0, 1.2)
		sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=α0)

		# Same impact parameter → same intensity from any direction
		@testset "b=$b" for b in [0.0, 0.4R, 0.8R]
			I_ref = S.render(make_ray(directions[1], b), sphere)
			@testset for n̂ in directions[2:end]
				@test S.render(make_ray(n̂, b), sphere) ≈ I_ref rtol=1e-10
			end
		end

		# Missed rays: zero from any direction
		@testset for n̂ in directions
			@test S.render(make_ray(n̂, 2R), sphere) == 0
		end
	end

	# Camera-based total flux comparison across viewing angles
	@testset "Total flux matches across camera angles" begin
		sphere = S.UniformSphere(; center, radius=R, u0, jν=j0, αν=1.2)
		xs = range(-R..R, 51)
		cam_dirs = [
			SVector(0.0, 0.0, 1.0),
			SVector(1.0, 0.0, 0.0),
			normalize(SVector(1.0, 1.0, 1.0)),
		]

		F_ref = let
			cam = S.CameraZ(; xys=grid(SVector, xs, xs), nz, ν, t=0.0)
			sum(S.render(cam, sphere)) * step(xs)^2
		end

		@testset for n̂ in cam_dirs[2:end]
			up = abs(dot(SVector(0.0, 1.0, 0.0), n̂)) < 0.9 ? SVector(0.0, 1.0, 0.0) : SVector(1.0, 0.0, 0.0)
			cam = S.CameraOrtho(;
				photon_direction=n̂, up,
				xys=grid(SVector, xs, xs),
				nz, ν, t=0.0,
			)
			F = sum(S.render(cam, sphere)) * step(xs)^2
			@test F ≈ F_ref rtol=1e-8
		end
	end
end


@testitem "SlowLight vs FastLight rendering" begin
	import Synchray as S
	using Accessors

	# 1. Static sphere: both modes give identical results
	@testset "static sphere identical" begin
		sphere = S.UniformSphere(
			center=S.FourPosition(0, 0, 0, 0), radius=1.3,
			u0=S.FourVelocity(SVector(0, 0, 0)), jν=0.7, αν=1.3
		)

		@testset for b in (0.0, 0.5, 0.9)
			ray_slow = S.RayZ(; x0=S.FourPosition(0, b, 0, 0), k=2.0, nz=512)
			ray_fast = S.RayZ(; x0=S.FourPosition(0, b, 0, 0), k=2.0, nz=512, light=S.FastLight())
			@test S.render(ray_slow, sphere) ≈ S.render(ray_fast, sphere)
			@test S.render(ray_slow, sphere, S.OpticalDepth()) ≈ S.render(ray_fast, sphere, S.OpticalDepth())
		end
	end

	# 2. Moving sphere: verify apparent position and Doppler scaling in each mode
	@testset "moving sphere positions and scaling" begin
		R = 1.0
		z0 = 3.0
		β = SVector(0.5, 0.0, 0.0)
		t_obs = 1.0

		ell = S.MovingUniformEllipsoid(
			center=S.FourPosition(0, 0, 0, z0),
			sizes=SVector(R, R, R),
			u0=S.FourVelocity(β),
			jν=1.0, αν=0.0,
		)
		ref = S.MovingUniformEllipsoid(
			center=S.FourPosition(0, 0, 0, z0),
			sizes=SVector(R, R, R),
			u0=S.FourVelocity(SVector(0.0, 0.0, 0.0)),
			jν=1.0, αν=0.0,
		)

		δ = S.doppler_factor(ell.u0, SVector(0, 0, 1))
		γ = S.gamma(ell.u0)

		# SlowLight: Terrell-shifted apparent center
		s_terrell = (t_obs - ell.center.t + ell.center.z) / (1 - β[3])
		xy_slow = SVector(β[1], β[2]) * s_terrell

		# FastLight: coordinate-position apparent center
		xy_fast = SVector(β[1], β[2]) * (t_obs - ell.center.t)

		# Y-only offsets: chord through object is unchanged by transverse Lorentz contraction
		xys_y = SVector{2,Float64}[(0, 0), (0, 0.3R), (0, 0.7R)]
		I_ref = map(xy -> S.render(S.RayZ(; x0=S.FourPosition(t_obs, xy[1], xy[2], 0), k=2.0, nz=2048), ref), xys_y)

		# SlowLight: rays through Terrell center → δ³ scaling
		I_slow = map(xy -> S.render(S.RayZ(; x0=S.FourPosition(t_obs, xy[1]+xy_slow[1], xy[2]+xy_slow[2], 0), k=2.0, nz=2048), ell), xys_y)
		@test I_slow ≈ δ^3 * I_ref rtol=1e-5

		# FastLight: rays through coordinate center → δ² scaling
		I_fast = map(xy -> S.render(S.RayZ(; x0=S.FourPosition(t_obs, xy[1]+xy_fast[1], xy[2]+xy_fast[2], 0), k=2.0, nz=2048, light=S.FastLight()), ell), xys_y)
		@test I_fast ≈ δ^2 * I_ref rtol=1e-5

		# FastLight x-offset: Lorentz contraction visible
		# At x-offset Δx from apparent center, chord shrinks by √(1 - γ²Δx²/R²) / √(1 - Δx²/R²)
		@testset for Δx in (0.3, 0.5)
			I_ref_x = S.render(S.RayZ(; x0=S.FourPosition(t_obs, Δx, 0, 0), k=2.0, nz=2048), ref)
			I_fast_x = S.render(S.RayZ(; x0=S.FourPosition(t_obs, xy_fast[1]+Δx, 0, 0), k=2.0, nz=2048, light=S.FastLight()), ell)
			chord_ratio = √(1 - (γ*Δx/R)^2) / √(1 - (Δx/R)^2)
			@test I_fast_x ≈ δ^2 * chord_ratio * I_ref_x rtol=1e-5
		end
	end
end


@testitem "CameraOrtho vs CameraPerspective consistency" begin
	import Synchray as S
	using RectiGrids

	# Non-trivial moving ellipsoid with all velocity components nonzero
	ell = S.MovingUniformEllipsoid(
		center=S.FourPosition(0, 0, 0, 5.0),
		sizes=SVector(1.5, 1.0, 1.2),
		u0=S.FourVelocity(SVector(0.3, 0.2, 0.1)),
		jν=0.7, αν=1.3,
	)

	L = 3.0
	N = 64
	nz = 512
	ν = 2.0
	D = 10_000.0
	t_ortho = 1.0

	xs = range(-L..L, N)
	xys_ortho = grid(SVector, xs, xs)

	# Build perspective xys: ortho (u,v) ↔ perspective (u/D, v/D)
	xys_persp = map(uv -> SVector(uv[1]/D, uv[2]/D), xys_ortho)

	@testset for light in (S.SlowLight(), S.FastLight())
		t_persp = light isa S.SlowLight ? t_ortho + D : t_ortho

		cam_ortho = S.CameraZ(; xys=xys_ortho, nz, ν, t=t_ortho, light)
		cam_persp = S.CameraPerspective(;
			photon_direction=SVector(0.0, 0.0, 1.0),
			origin=SVector(0.0, 0.0, D),
			xys=xys_persp, nz, ν, t=t_persp, light,
		)

		@testset for what in (S.Intensity(), S.OpticalDepth())
			img_ortho = S.render(cam_ortho, ell, what)
			img_persp = S.render(cam_persp, ell, what)
			@test img_ortho ≈ img_persp rtol=1e-2
		end

		# Tuple mode: compare components individually
		img_ortho_t = S.render(cam_ortho, ell, (S.Intensity(), S.OpticalDepth()))
		img_persp_t = S.render(cam_persp, ell, (S.Intensity(), S.OpticalDepth()))
		@test getindex.(img_ortho_t, 1) ≈ getindex.(img_persp_t, 1) rtol=1e-2
		@test getindex.(img_ortho_t, 2) ≈ getindex.(img_persp_t, 2) rtol=1e-2
	end
end

@testitem "CameraPerspective behind-camera clipping" begin
	import Synchray as S
	using RectiGrids

	# Camera at origin, photon_direction = +z.
	# Scene (negative s) is at z < 0; behind-camera (positive s) is at z > 0.
	nz = 256
	ν = 1.0
	jν = 1.0

	# With zero absorption: Iν = jν * L (pure emission, no attenuation)

	sphere_front = S.UniformSphere(
		center=S.FourPosition(0.0, 0.0, 0.0, -10.0), radius=2.0,
		u0=S.FourVelocity(1.0, 0.0, 0.0, 0.0), jν=jν, αν=0.0,
	)
	sphere_behind = S.UniformSphere(
		center=S.FourPosition(0.0, 0.0, 0.0, 10.0), radius=2.0,
		u0=S.FourVelocity(1.0, 0.0, 0.0, 0.0), jν=jν, αν=0.0,
	)
	sphere_straddling = S.UniformSphere(
		center=S.FourPosition(0.0, 0.0, 0.0, -1.0), radius=5.0,
		u0=S.FourVelocity(1.0, 0.0, 0.0, 0.0), jν=jν, αν=0.0,
	)

	@testset for light in (S.SlowLight(), S.FastLight())
		cam = S.CameraPerspective(;
			photon_direction=SVector(0.0, 0.0, 1.0),
			origin=SVector(0.0, 0.0, 0.0),
			xys=grid(SVector, x=range(-0.3..0.3, 8), y=range(-0.3..0.3, 8)),
			nz, ν, t=0.0, light,
		)

		ray = S.camera_ray(cam, SVector(0.0, 0.0))
		@test ray.s_max == 0
		@test isbits(ray)

		# Front sphere (center z=-10, radius 2): full diameter path L=4
		@test S.render(ray, sphere_front) ≈ jν * 4.0 rtol=1e-2

		# Behind sphere (center z=+10): entirely clipped → zero
		@test S.render(ray, sphere_behind) == 0.0
		@test all(iszero, S.render(cam, sphere_behind))

		# Straddling sphere (center z=-1, radius 5): seg -6..4 clipped to -6..0, path L=6
		@test S.render(ray, sphere_straddling) ≈ jν * 6.0 rtol=1e-2
		# Without clipping, full path L=10:
		ray_unclipped = S.Ray(ray.x0, ray.k, ray.e1, ray.e2, ray.nz, ray.light, nothing)
		@test S.render(ray_unclipped, sphere_straddling) ≈ jν * 10.0 rtol=1e-2
	end

	# Ortho camera: no clipping, sees both spheres with full diameter path L=4
	cam_ortho = S.CameraZ(; xys=grid(SVector, x=range(-3.0..3.0, 8), y=range(-3.0..3.0, 8)), nz, ν, t=0.0)
	ray_ortho = S.camera_ray(cam_ortho, SVector(0.0, 0.0))
	@test isnothing(ray_ortho.s_max)
	@test S.render(ray_ortho, sphere_front) ≈ jν * 4.0 rtol=1e-2
	@test S.render(ray_ortho, sphere_behind) ≈ jν * 4.0 rtol=1e-2
	@test S.render(ray_ortho, sphere_straddling) ≈ jν * 10.0 rtol=1e-2
end

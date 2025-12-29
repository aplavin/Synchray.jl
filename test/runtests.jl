using TestItems
using TestItemRunner
@run_package_tests


@testitem "Minkowski + Doppler conventions" begin
	import Synchray as S

	@testset "arithmetics" begin
		a = S.FourPosition(1, 2, 3, 4)
		b = S.FourPosition(0.5, -1.0, 2.0, 3.0)
		@test a + b === S.FourPosition(1.5, 1.0, 5.0, 7.0)
		@test a - b === S.FourPosition(0.5, 3.0, 1.0, 1.0)
		@test 2a === S.FourPosition(2, 4, 6, 8)
	end

	@testset "minkowski_dot basics" begin
		a = S.FourPosition(2.0, 1.0, 3.0, -4.0)
		b = S.FourPosition(-1.0, 2.0, 0.5, 7.0)
		@test S.minkowski_dot(a, b) ≈ S.minkowski_dot(b, a)
		@test S.minkowski_dot(a, a) ≈ -(a.t^2) + a.x^2 + a.y^2 + a.z^2
	end

	@testset "FourVelocity normalization" begin
		β = SVector(0.3, -0.4, 0.1)
		u = S.FourVelocity(β)
		@test S.minkowski_dot(u, u) ≈ -1 atol=1e-12
		@test u.t ≈ inv(sqrt(1 - dot(β, β)))
		@test S.beta(u) ≈ β atol=1e-12
		@test S.gamma(u) ≈ u.t atol=1e-12
		@test S.gamma(β) ≈ u.t atol=1e-12
	end

	@testset "photon_k is null" begin
		ν = 2.5
		n = normalize(SVector(0.2, -0.3, 0.7))
		k = S.photon_k(ν, n)
		@test S.minkowski_dot(k, k) ≈ 0 atol=1e-12
	end

	@testset "Doppler boosting δ convention" begin
		# Convention locked in by this test:
		# doppler_factor(u, n) == δ = ν_obs / ν_comoving
		# For fast motion toward the observer along the photon direction (β ⋅ n > 0): δ >> 1.
		β = 0.99
		u_toward = S.FourVelocity(SVector(0.0, 0.0, β))
		n = SVector(0.0, 0.0, 1.0)

		δ_toward = S.doppler_factor(u_toward, n)
		@test δ_toward > 1
		@test δ_toward > 10

		u_away = S.FourVelocity(SVector(0.0, 0.0, -β))
		δ_away = S.doppler_factor(u_away, n)
		@test δ_away < 1
	end
end


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
				ℓ = 2 * sqrt(R^2 - b^2)
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
				ℓ = t > 0 ? 2 * sqrt(t / A) : 0.0
				(δ^3) * ℓ
			end, xys)
			@test count(>(0), Iexp) == 4

			@test Icur ≈ Iexp rtol=3e-3
		end
	end
end


@testitem "Synchrotron slab scalings" begin
	import Synchray as S
	using Accessors

	L = 2
	ν = 2.0
	model = S.PowerLawElectrons(; p=3.2, Cj=0.7, Ca=0.0)
	ne0 = 1.3
	B0 = 0.9
	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=ν, nz=4_000)

	cases = (
		S.FourVelocity(SVector(0, 0, 0)),
		S.FourVelocity(SVector(0, 0, 0.3)),
		S.FourVelocity(SVector(0.3, 0, 0.5)),
	)

	@testset for u0 in cases
		δ = S.doppler_factor(u0, SVector(0, 0, 1))
		ν′ = ν / δ

		slab = S.UniformSynchrotronSlab(; z=0..L, u0, ne0, B0, model)

		@testset "optically thin" begin
			I_num = S.render(ray, slab)
			@test I_num > 0
			@test S.render(ray, slab, S.OpticalDepth()) ≈ 0 atol=1e-12

			(j0, _) = S._synchrotron_coeffs(model, ne0, B0, ν′)
			I_exact = j0 * L * (δ^2)
			@test I_num ≈ I_exact rtol=2e-3

			@test S.render(ray, @set slab.ne0 *= 2) ≈ 2I_num rtol=2e-3
			@test S.render(ray, @set slab.B0 *= 2) ≈ I_num * (2^((model.p + 1) / 2)) rtol=2e-3
			@test S.render((@set S.photon_frequency(ray) *= 2), slab) ≈ I_num * (2^(-(model.p - 1) / 2)) rtol=2e-3

			αobs = S.render(ray, slab, S.SpectralIndex())
			@test αobs ≈ (-(model.p - 1) / 2) atol=2e-4
		end

		@testset "very optically thick" begin
			model_thick = @set model.Ca = 5e6
			slab_thick = @set slab.model = model_thick
			(j_thick, α_thick) = S._synchrotron_coeffs(model_thick, ne0, B0, ν′)

			I_num_thick = S.render(ray, slab_thick)
			I_exact_thick = (j_thick / α_thick) * (δ^3)
			@test I_num_thick ≈ I_exact_thick rtol=2e-3
			@test S.render(ray, slab_thick, S.OpticalDepth()) ≈ (α_thick * (L / δ)) rtol=2e-3
			@test S.render(ray, slab_thick, S.SpectralIndex()) ≈ 5/2 atol=2e-4
			@test S.render((@set S.photon_frequency(ray) *= 2), slab_thick) ≈ I_num_thick * (2^(5/2)) rtol=2e-3

			@test S.render(ray, @set slab_thick.ne0 *= 2) ≈ I_num_thick rtol=1e-3
		end
	end
end


@testitem "Conical BK jet geometry + scalings" begin
	import Synchray as S
	using Accessors

	jet = S.ConicalBKJet(; 
		axis=SVector(0, 0, 1),
		φj=0.05,
		s=0.1..5,
		s0=1,
		ne0=2,
		B0=3,
		speed_profile=(η -> 0),
		model=S.PowerLawElectrons(; p=2.5, Cj=1, Ca=1),
	)

	# Scalings along the axis at s=s0 and s=2s0.
	x4_1 = S.FourPosition(0, 0, 0, 1)
	x4_2 = S.FourPosition(0, 0, 0, 2)
	@test S.electron_density(jet, x4_1) ≈ 2 atol=1e-12
	@test S.magnetic_field_strength(jet, x4_1) ≈ 3 atol=1e-12
	@test S.electron_density(jet, x4_2) ≈ 2 * (2^(-2)) atol=1e-12
	@test S.magnetic_field_strength(jet, x4_2) ≈ 3 * (2^(-1)) atol=1e-12

	# Outside cone should yield zero microphysics.
	x4_out = S.FourPosition(0, 3 * 2 * tan(jet.φj), 0, 2)
	@test S.electron_density(jet, x4_out) == 0
	@test S.magnetic_field_strength(jet, x4_out) == 0

	ray = S.RayZ(; x0=S.FourPosition(0, 0, 0, 0), k=2, nz=1024)
	# On-axis ray should see the full truncation segment for axis=ẑ.
	@test S.z_interval(jet, ray) == jet.s

	# A ray outside the projected cone should miss entirely.
	xmiss = 1.1 * maximum(jet.s) * tan(jet.φj)
	miss_ray = @set ray.x0.x = xmiss
	@test S.render(miss_ray, jet) == 0
	@test S.render(miss_ray, jet, S.OpticalDepth()) == 0
	@test S.render(miss_ray, jet, S.SpectralIndex()) |> isnan
end


# @testitem "_" begin
#     import Aqua
#     Aqua.test_all(Synchray)

#     import CompatHelperLocal as CHL
#     CHL.@check()
# end

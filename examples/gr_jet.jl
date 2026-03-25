include("helpers.jl")
using Krang
using StructArrays


function gr_jet_image(; log_path, suffix="")
	log_section!(log_path, "gr_jet_image")

	jet = let
		axis = SVector(0.0, 0.0, 1.0)
		S.EmissionRegion(;
			geometry = S.Geometries.Conical(; axis, φj=15u"°", z=5.0..200.0),
			velocity = S.VelocitySpec(S.Directions.Axial(), S.gamma, S.Profiles.Constant(20.0)),
			emission = S.SynchrotronEmission(
				ne = S.Profiles.Axial(S.PowerLaw(-2; val0=100.0, s0=1.0)),
				B = S.BFieldSpec(
					S.Profiles.Axial(S.PowerLaw(-1; val0=500.0, s0=1.0)),
					S.Directions.Scalar(),
					b -> S.FullyTangled(b),
				),
				electrons = S.IsotropicPowerLawElectrons(; p=2.5, Cj=1.0, Ca=1.0),
			),
		) |> S.prepare_for_computations
	end

	spin = 0.5
	camera_extent = 60.0
	npx = 256
	nz = 256
	s_range_3d = -200.0 .. 200.0

	views = [
		(name="5", θ=5u"°"),
		(name="20", θ=20u"°"),
		(name="45", θ=45u"°"),
		(name="90", θ=90u"°"),
		(name="160", θ=160u"°"),
		(name="175", θ=175u"°"),
	]

	# render all views first for shared colorscale
	results = map(views) do v
		n = SVector(sin(v.θ), 0.0, cos(v.θ))
		cam_flat = S.CameraOrtho(;
			photon_direction = n,
			xys = grid(SVector,
				x = range(-camera_extent, camera_extent, length=npx) |> LinRange,
				y = range(-camera_extent, camera_extent, length=npx) |> LinRange),
			nz, ν=1e11, t=0.0,
		)
		lens = S.GRLens(; spin)
		cam_gr = S.CameraGR(; camera=cam_flat, lens)
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam_gr, jet)
				label = f"GR jet, θ={v.name}°, {npx}px × {nz}"
				@time label S.render(cam_gr, jet)
			end
		end
		(; img, cam_flat, lens)
	end

	global_max = maximum(r -> maximum(r.img), results)

	fig = Figure()

	for (col, (v, r)) in enumerate(zip(views, results))
		# top row: rendered GR image
		ax = Axis(fig[1, col];
			title=f"θ = {v.name}°, a = {spin}",
			aspect=DataAspect(), width=200, height=200)
		plt = image!(ax, r.img; colormap=:inferno, colorscale=SymLog(1e-10 * global_max))
		hidespines!(ax)
		hidedecorations!(ax)
		col == length(views) && Colorbar(fig[1, col+1], plt; tickformat=EngTicks())

		# bottom row: 3D ray paths
		ext3d = 100.0
		ax3d = Axis3(fig[2, col];
			limits=((-ext3d, ext3d), (-ext3d, ext3d), (-ext3d, ext3d)),
			aspect=:data,
			title=f"θ = {v.name}°")

		# BH marker
		bh_local = S.rotate_lab_to_local(jet.geometry, r.lens.bh_position)
		meshscatter!(ax3d, [Point3f(bh_local...)]; markersize=2, color=:black)

		# sample rays: 3x3 grid across the FOV
		sample_uvs = [SVector(u, vv) for u in LinRange(-30.0, 30.0, 3) for vv in LinRange(-30.0, 30.0, 3)]
		for uv in sample_uvs
			ray_in = S.camera_ray(r.cam_flat, uv)
			raygr = S.RayGR2(ray_in, r.lens)
			# straight ray (dashed)
			flat_pts = S.ray_in_local_coords(ray_in, jet.geometry; s_range=s_range_3d)
			linesegments!(ax3d, [Point3f(p...) for p in flat_pts]; linestyle=:dash, linewidth=1.5, color=Cycled(1))
			# GR-deflected ray (solid, colored by segment)
			gr_pts = S.ray_in_local_coords(raygr, jet.geometry; s_range=s_range_3d)
			gr_colors = Makie.wong_colors()[[1,1,2,2]][1:length(gr_pts)]
			linesegments!(ax3d, [Point3f(p...) for p in gr_pts]; color=gr_colors, linewidth=1.5)
		end
	end

	resize_to_layout!(fig)
	save(joinpath(outdir, "gr_jet$(suffix).png"), fig)
	fig
end


function main(; log_path, suffix="")
	gr_jet_image(; log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
	logs = init_script(@__FILE__)
	main(; log_path=logs.log_cpu)
end

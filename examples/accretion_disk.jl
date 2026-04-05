include("helpers.jl")
using Krang
using FastChebInterp
using OhMyThreads: tmap

_tmap(f, X::KeyedArray) = @set AxisKeys.keyless_unname(X) = tmap(f, X)


function accretion_disk_image(; log_path, suffix="")
	log_section!(log_path, "accretion_disk_image")

	spin = 0.94

	disk = let
		S.EmissionRegion(;
			geometry = S.Geometries.PowerLawDisk(r_range=3.0..20.0, h_ref=0.5, r_ref=1.0, a=0),
			velocity = S.KeplerianVelocity(spin=spin),
			emission = S.FixedEmission(S=1.0, α=0.5),
		) |> S.prepare_for_computations
	end

	npx = 256
	camera_extent = 25.0

	views = [
		(name="15", θ=15u"°"),
		(name="75", θ=75u"°"),
		(name="85", θ=85u"°"),
	]

	results = map(views) do v
		n = SVector(sin(v.θ), 0.0, cos(v.θ))
		cam_flat = S.CameraOrtho(;
			photon_direction = n,
			xys = grid(SVector,
				x = range(-camera_extent, camera_extent, length=npx) |> LinRange,
				y = range(-camera_extent, camera_extent, length=npx) |> LinRange),
			nz = 8, ν = 1.0, t = 0.0,
			mapfunc = _tmap,
		)
		cam = S.CameraKerrGR(; camera=cam_flat, metric_spin=spin, nτ=800, τ_range=(0.01, 0.99))
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam, disk)
				label = f"Accretion disk, θ={v.name}°, {npx}px"
				@time label S.render(cam, disk)
			end
		end
		(; img)
	end

	global_max = maximum(r -> maximum(r.img), results)

	fig = Figure(size=(1600, 500))

	for (col, (v, r)) in enumerate(zip(views, results))
		ax = Axis(fig[1, col];
			title=f"θ = {v.name}°, a = {spin}",
			aspect=DataAspect(), width=450, height=450)
		plt = image!(ax, r.img; colormap=:inferno)
		hidespines!(ax)
		hidedecorations!(ax)
		col == length(views) && Colorbar(fig[1, col+1], plt; tickformat=EngTicks())
	end

	resize_to_layout!(fig)
	save(joinpath(outdir, "accretion_disk$(suffix).png"), fig)
	fig
end


function accretion_disk_image_cached(; log_path, suffix="")
	log_section!(log_path, "accretion_disk_image_cached")

	spins = [0.1, 0.98]

	npx = 256
	camera_extent = 10.0

	views = [
		(name="15", θ=15u"°"),
		(name="75", θ=75u"°"),
		(name="85", θ=85u"°"),
	]

	make_disk(spin) = S.EmissionRegion(;
		geometry = S.Geometries.PowerLawDisk(r_range=2.0..7.0, h_ref=0.5, r_ref=1.0, a=0),
		velocity = S.KeplerianVelocity(spin=spin),
		emission = S.FixedEmission(S=1.0, α=0.5),
	) |> S.prepare_for_computations

	make_cam_flat(v) = let n = SVector(sin(v.θ), 0.0, cos(v.θ))
		S.CameraOrtho(;
			photon_direction = n,
			xys = grid(SVector,
				x = range(-camera_extent, camera_extent, length=npx) |> LinRange,
				y = range(-camera_extent, camera_extent, length=npx) |> LinRange),
			nz = 8, ν = 1.0, t = 0.0,
			mapfunc = _tmap,
		)
	end

	# Row 1: flat space (use first spin for Keplerian velocity)
	disk_flat = make_disk(first(spins))
	flat_results = map(views) do v
		cam_flat = make_cam_flat(v)
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam_flat, disk_flat)
				label = f"Accretion disk flat, θ={v.name}°, {npx}px"
				@time label S.render(cam_flat, disk_flat)
			end
		end
		(; img)
	end

	# Rows 2+: GR cached at different spins
	gr_results = map(spins) do spin
		disk = make_disk(spin)
		map(views) do v
			cam_flat = make_cam_flat(v)
			cam = S.CameraKerrGRCached(; camera=cam_flat, metric_spin=spin, nτ=400, R=10.0)
			img = open(log_path, "a") do io
				redirect_stdout(io) do
					S.render(cam, disk)
					label = f"Accretion disk cached, a={spin}, θ={v.name}°, {npx}px"
					@time label S.render(cam, disk)
				end
			end
			(; img)
		end
	end

	fig = Figure(size=(1600, 1500))

	for (col, (v, r)) in enumerate(zip(views, flat_results))
		ax = Axis(fig[1, col];
			title=f"θ = {v.name}°, flat",
			aspect=DataAspect(), width=450, height=450)
		plt = image!(ax, r.img; colormap=:inferno)
		hidespines!(ax)
		hidedecorations!(ax)
		col == length(views) && Colorbar(fig[1, col+1], plt; tickformat=EngTicks())
	end

	for (row, (spin, results)) in enumerate(zip(spins, gr_results))
		for (col, (v, r)) in enumerate(zip(views, results))
			ax = Axis(fig[row+1, col];
				title=f"θ = {v.name}°, a = {spin}",
				aspect=DataAspect(), width=450, height=450)
			plt = image!(ax, r.img; colormap=:inferno)
			hidespines!(ax)
			hidedecorations!(ax)
			col == length(views) && Colorbar(fig[row+1, col+1], plt; tickformat=EngTicks())
		end
	end

	resize_to_layout!(fig)
	save(joinpath(outdir, "accretion_disk_cached$(suffix).png"), fig)
	fig
end


function conical_jet_image(; log_path, suffix="")
	log_section!(log_path, "conical_jet_image")

	spin = 0.94

	jet = let
		S.EmissionRegion(;
			geometry = S.Geometries.Conical(axis=SVector(0, 0, 1), φj=deg2rad(60), z=1.0..20.0),
			velocity = S.VelocitySpec(S.Directions.Axial(), S.Profiles.Constant(10.0)),
			emission = S.FixedEmission(S=1.0, α=0.005),
		) |> S.prepare_for_computations
	end

	npx = 256
	camera_extent = 25.0

	views = [
		(name="15", θ=15u"°"),
		(name="75", θ=75u"°"),
		(name="85", θ=85u"°"),
	]

	make_cam_flat(v) = let n = SVector(sin(v.θ), 0.0, cos(v.θ))
		S.CameraOrtho(;
			photon_direction = n,
			xys = grid(SVector,
				x = range(-camera_extent, camera_extent, length=npx) |> LinRange,
				y = range(-camera_extent, camera_extent, length=npx) |> LinRange),
			nz = 8, ν = 1.0, t = 0.0,
			mapfunc = _tmap,
		)
	end

	# Row 1: flat space
	flat_results = map(views) do v
		cam_flat = make_cam_flat(v)
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam_flat, jet)
				label = f"Conical jet flat, θ={v.name}°, {npx}px"
				@time label S.render(cam_flat, jet)
			end
		end
		(; img)
	end

	# Row 2: GR cached
	gr_results = map(views) do v
		cam_flat = make_cam_flat(v)
		cam = S.CameraKerrGRCached(; camera=cam_flat, metric_spin=spin, nτ=400, R=21.0)
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam, jet)
				label = f"Conical jet GR, θ={v.name}°, {npx}px"
				@time label S.render(cam, jet)
			end
		end
		(; img)
	end

	fig = Figure(size=(1600, 1000))

	for (col, (v, r)) in enumerate(zip(views, flat_results))
		ax = Axis(fig[1, col];
			title=f"θ = {v.name}°, flat",
			aspect=DataAspect(), width=450, height=450)
		plt = image!(ax, r.img; colormap=:inferno)
		hidespines!(ax)
		hidedecorations!(ax)
		col == length(views) && Colorbar(fig[1, col+1], plt; tickformat=EngTicks())
	end

	for (col, (v, r)) in enumerate(zip(views, gr_results))
		ax = Axis(fig[2, col];
			title=f"θ = {v.name}°, a = {spin}",
			aspect=DataAspect(), width=450, height=450)
		plt = image!(ax, r.img; colormap=:inferno)
		hidespines!(ax)
		hidedecorations!(ax)
		col == length(views) && Colorbar(fig[2, col+1], plt; tickformat=EngTicks())
	end

	resize_to_layout!(fig)
	save(joinpath(outdir, "conical_jet$(suffix).png"), fig)
	fig
end


function main(; log_path, suffix="")
	accretion_disk_image(; log_path, suffix)
	accretion_disk_image_cached(; log_path, suffix)
	conical_jet_image(; log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
	logs = init_script(@__FILE__)
	main(; log_path=logs.log_cpu)
end

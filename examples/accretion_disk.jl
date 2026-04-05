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

	npx = 512
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
			nz = 1, ν = 1.0, t = 0.0,
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

	spin = 0.94

	disk = let
		S.EmissionRegion(;
			geometry = S.Geometries.PowerLawDisk(r_range=3.0..20.0, h_ref=0.5, r_ref=1.0, a=0),
			velocity = S.KeplerianVelocity(spin=spin),
			emission = S.FixedEmission(S=1.0, α=0.5),
		) |> S.prepare_for_computations
	end

	npx = 512
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
			nz = 1, ν = 1.0, t = 0.0,
			mapfunc = _tmap,
		)
		cam = S.CameraKerrGRCached(; camera=cam_flat, metric_spin=spin, nτ=800, τ_range=(0.01, 0.99), R=camera_extent)
		img = open(log_path, "a") do io
			redirect_stdout(io) do
				S.render(cam, disk)
				label = f"Accretion disk cached, θ={v.name}°, {npx}px"
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
	save(joinpath(outdir, "accretion_disk_cached$(suffix).png"), fig)
	fig
end


function main(; log_path, suffix="")
	accretion_disk_image(; log_path, suffix)
	accretion_disk_image_cached(; log_path, suffix)
end

if abspath(PROGRAM_FILE) == @__FILE__
	logs = init_script(@__FILE__)
	main(; log_path=logs.log_cpu)
end

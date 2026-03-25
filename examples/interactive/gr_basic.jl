using MakieExtra; import GLMakie
using RectiGrids, AxisKeysExtra
using AccessibleModels; using AccessibleModels: P
using Rotations
using DistributionsExtra
using PyFormattedStrings
using DataManipulation
using StructArrays

using Synchray
import Synchray as S
using Krang

MakieExtra.show_gl_icon_in_dock()

# ── Fixed jet parameters ──
jet = let
	axis = SVector(0.0, 0.0, 1.0)
	S.EmissionRegion(;
		geometry = S.Geometries.Conical(; axis, φj=deg2rad(15), z=5.0..200.0),
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

# ── UI ──
fig = Figure(size=(1200, 1000))

params, = SliderGrid(GridLayout(fig[1, 1], tellheight=false)[1,1], AccessibleModel((;
	viewing_ang = P(1..179, 20)u"°",
	spin = P(LogUniform(0.01..0.998), 0.5),
	npix = P(discreterange(log, 32..300, length=30), 80),
	nz = P(discreterange(log, 10..500, length=30), 100),
	orthogonal = P([false, true], true),
), AccessibleModels.Auto()))

# ── Camera setup ──
camera_label(use_ortho) = use_ortho ? "Ortho" : "Perspective"

camera_distance = 1000.0
camera_extent = 60.0

camera_ortho = @lift let
	θ = $params.viewing_ang
	n = RotY(θ) * SVector(0.0, 0.0, 1.0)  # photon direction (from BH toward observer)
	S.CameraOrtho(;
		photon_direction = n,
		xys = grid(SVector,
			x = range(-camera_extent, camera_extent, length=$params.npix) |> LinRange,
			y = range(-camera_extent, camera_extent, length=$params.npix) |> LinRange),
		nz = $params.nz,
		ν = 1e11,
		t = 0.0,
	)
end

camera_perspective = @lift let
	θ = $params.viewing_ang
	n = RotY(θ) * SVector(0.0, 0.0, 1.0)  # photon direction (from BH toward observer)
	S.CameraPerspective(;
		photon_direction = n,
		origin = camera_distance * n,  # far from BH, looking toward it
		xys = grid(SVector,
			x = range(-camera_extent / camera_distance, camera_extent / camera_distance, length=$params.npix) |> LinRange,
			y = range(-camera_extent / camera_distance, camera_extent / camera_distance, length=$params.npix) |> LinRange),
		nz = $params.nz,
		ν = 1e11,
		t = 0.0,
	)
end

lens = @lift S.GRLens(; spin=$params.spin)

flat_camera = @lift($params.orthogonal ? $camera_ortho : $camera_perspective)::Any

gr_camera = @lift S.CameraGR(; camera=$flat_camera, lens=$lens)

# ── Render ──
img_flat = @lift S.render($flat_camera, jet)
img_gr = @lift S.render($gr_camera, jet)

# ── Plot ──
colorscale = @lift SymLog(1e-3 * max(maximum($img_flat), maximum($img_gr)))

ax_flat = Axis(fig[1, 2];
	title = @lift(f"{camera_label($params.orthogonal)} jet  θ={$params.viewing_ang:.0f}"),
	xlabel = "x [rg]", ylabel = "y [rg]",
	aspect = DataAspect())
plt_flat = image!(ax_flat, img_flat; colorscale, colormap = :inferno)

ax_gr = Axis(fig[1, 3];
	title = @lift(f"GR  a={$params.spin:.2f}  θ={$params.viewing_ang:.0f}"),
	xlabel = "x [rg]", ylabel = "y [rg]",
	aspect = DataAspect())
plt_gr = image!(ax_gr, img_gr; colorscale, colormap = :inferno)

Colorbar(fig[1, 4], plt_gr, tickformat=EngTicks())

# ── Mouse selection ──
uv_sel = Observable(SVector(30.0, 0.0))

on(mouse_position_obs(ax_flat; key=Mouse.left)) do pos
	uv_sel[] = SVector(pos...)
end
on(mouse_position_obs(ax_gr; key=Mouse.left)) do pos
	uv_sel[] = SVector(pos...)
end

# Mark selected pixel on both images
scatter!(ax_flat, @lift(Point2f($uv_sel...)); color=:transparent, markersize=15, marker=:circle, strokewidth=1, strokecolor=:cyan)
scatter!(ax_gr, @lift(Point2f($uv_sel...)); color=:transparent, markersize=15, marker=:circle, strokewidth=1, strokecolor=:cyan)

# # ── 3D ray visualization ──
ray_gr_sel = @lift let
	uv = $uv_sel
	cam = $flat_camera
	l = $lens
	ray_in = S.camera_ray(cam, uv)
	S.RayGR2(ray_in, l)
end::Any

s_range_3d = -500.0 .. 500.0

ray_path_3d = @lift let
	ray = $ray_gr_sel
	geom = jet.geometry
	pts = S.ray_in_local_coords(ray, geom; s_range=s_range_3d)
	colors = Makie.wong_colors()[[1,1,2,2]][1:length(pts)]
	StructArray(;pts, colors)
end

ax3d = Axis3(fig[2, 2:3];
	limits=((-20, 20), (-20, 20), (-20, 20)),
	xlabel="x", ylabel="y", zlabel="z (jet axis)",
	aspect=:data,
	title="GR ray path (local coords)")

# BH sphere
bh_local = @lift let
	S.rotate_lab_to_local(jet.geometry, $lens.bh_position)
end
meshscatter!(ax3d, @lift(Point3f($bh_local...)); markersize=1, color=:black)

# GR ray path
linesegments!(ax3d, (@lift FPlot($ray_path_3d, first, color=last)); linewidth=2)


display(fig)

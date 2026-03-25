using MakieExtra; import GLMakie
using RectiGrids, AxisKeysExtra
using AccessibleModels; using AccessibleModels: P
using Rotations
using DistributionsExtra
using PyFormattedStrings
using DataManipulation

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
fig = Figure(size=(1200, 700))

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

is_perspective = @lift !$params.orthogonal

plot_xs_flat = @lift let
	xs = collect(Float64, axiskeys($img_flat, 1))
	$is_perspective ? camera_distance .* xs : xs
end
plot_ys_flat = @lift let
	ys = collect(Float64, axiskeys($img_flat, 2))
	$is_perspective ? camera_distance .* ys : ys
end
plot_xs_gr = @lift let
	xs = collect(Float64, axiskeys($img_gr, 1))
	$is_perspective ? camera_distance .* xs : xs
end
plot_ys_gr = @lift let
	ys = collect(Float64, axiskeys($img_gr, 2))
	$is_perspective ? camera_distance .* ys : ys
end

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

display(fig)

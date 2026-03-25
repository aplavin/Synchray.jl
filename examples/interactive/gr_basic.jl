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
		geometry = S.Geometries.Conical(; axis, φj=deg2rad(15), z=20.0..200.0),
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
	gr = P([false, true], true),
), AccessibleModels.Auto()))

# ── Camera setup ──
camera_label(use_ortho, use_gr) = use_gr ? "GR" : use_ortho ? "Ortho" : "Perspective"

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

lens = @lift let
	S.GRLens(; spin=$params.spin)
end

flat_camera = @lift let
	$params.orthogonal ? $camera_ortho : $camera_perspective
end

defl_map = @lift let
	S.compute_deflection_map($lens, $flat_camera)
end

camera = @lift let
	if $params.gr
		S.CameraGR(; camera=$flat_camera, deflection=$defl_map, lens=$lens)
	else
		$flat_camera
	end
end::Any

# ── Render ──
img = @lift S.render($camera, jet)

uses_perspective(cam) =
	cam isa S.CameraPerspective || (cam isa S.CameraGR && cam.camera isa S.CameraPerspective)

plot_xs = @lift let
	xs = collect(Float64, axiskeys($img, 1))
	uses_perspective($camera) ? camera_distance .* xs : xs
end

plot_ys = @lift let
	ys = collect(Float64, axiskeys($img, 2))
	uses_perspective($camera) ? camera_distance .* ys : ys
end

# ── Plot ──
ax = Axis(fig[1, 2];
	title = @lift(f"{camera_label($params.orthogonal, $params.gr)} jet  θ={$params.viewing_ang:.0f}  a={$params.spin:.2f}"),
	xlabel = "x [rg]", ylabel = "y [rg]",
	aspect = DataAspect())

plt = image!(ax, img; colorscale=(@lift SymLog(1e-3*maximum($img))),colormap = :inferno)
Colorbar(fig[1, 3], plt, tickformat=EngTicks())

display(fig)

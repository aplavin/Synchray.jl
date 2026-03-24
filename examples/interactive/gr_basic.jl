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
		velocity = S.VelocitySpec(S.Directions.Axial(), S.gamma, S.Profiles.Constant(3.0)),
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
), AccessibleModels.Auto()))

# ── Camera: rotates around the jet ──
camera = @lift let
	θ = $params.viewing_ang
	n = RotY(θ) * SVector(0.0, 0.0, -1.0)  # look direction toward BH
	S.CameraOrtho(;
		look_direction = n,
		origin = -1000.0 * n,  # far from BH, looking toward it
		xys = grid(SVector,
			x = range(-60.0, 60.0, length=$params.npix) |> LinRange,
			y = range(-60.0, 60.0, length=$params.npix) |> LinRange),
		nz = $params.nz,
		ν = 1e11,
		t = 0.0,
	)
end

# ── Deflection map (only depends on viewing angle + spin) ──
defl_map = @lift let
	θ_obs = Float64($params.viewing_ang)
	# Krang needs θ_obs away from 0 and π
	θ_obs_clamped = clamp(θ_obs, 0.02, π - 0.02)
	S.compute_deflection_map($params.spin, θ_obs_clamped, $camera)
end

# ── Render ──
img = @lift let
	cam_gr = S.CameraGR(; camera=$camera, deflection=$defl_map)
	S.render(cam_gr, jet)
end

# ── Plot ──
ax = Axis(fig[1, 2];
	title = @lift(f"GR jet  θ={$params.viewing_ang:.0f}  a={$params.spin:.2f}"),
	xlabel = "x [rg]", ylabel = "y [rg]",
	aspect = DataAspect())

heatmap!(ax,
	@lift(collect(Float64, axiskeys($img, 1))),
	@lift(collect(Float64, axiskeys($img, 2))),
	@lift(log10.(max.(collect($img), 1e-30)));
	colormap = :inferno)

display(fig)

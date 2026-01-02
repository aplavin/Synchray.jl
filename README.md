# Synchray.jl

Special-relativistic raytracing and radiative transfer for synchrotron emission, aimed at studing astrophysical jets. Renders images with an orthographic camera.

---

Current content: SR, straight rays, synchrotron Stokes I, variety of built-in objects.

Designed for flexibility and extensibility, easy to define new objects, geometries, and emission models.

Fast, suitable for interactive use and parameter exploration.

## Install

```julia
] add https://github.com/JuliaAPlavin/Synchray.jl
```

## Minimal use (conical jet)

```julia
import Synchray as S

# Image-plane grid
xs = range(-1.0, 1.0; length=256)
ys = range(-1.0, 1.0; length=256)
cam = S.CameraZ(
    xys=S.grid(S.SVector, xs, ys),
    nz=200, ν=230.0, t=0.0)

jet = S.ConicalBKJet(
	axis = SVector(0.0, 0.0, 1.0),
	φj = 0.03,
	s = 0.1..10.0,
	s0 = 1.0,
	ne0 = 1.0,
	B0 = 1.0,
	speed_profile = η -> (S.beta, 0.99),
)

Iν = S.render(cam, jet, S.Intensity())
# Iν is a 256×256 array of per-pixel specific intensities
```

`Unitful` inputs/outputs are supported: use `S.withunits(S.ConicalBKJet; ...)`, `S.withunits(S.CameraZ; ...)`, and `S.withunits(S.render, cam, jet, S.Intensity())`.

# Synchray.jl

Special-relativistic raytracing and radiative transfer for synchrotron emission, aimed at studing astrophysical jets.


https://github.com/user-attachments/assets/7ef2649c-6297-4fa8-983c-9ea201995712


*Code for this example: https://plav.in/Synchray/viewing_angle (notebook).*

Current content: volumetric radiative transfer, SR, synchrotron Stokes I and linear polarization, autodiff, variety of built-in objects.

Designed for flexibility and extensibility, easy to define new objects, geometries, and emission models.

Fast – HD images in 60 FPS on a laptop – suitable for interactive use and parameter exploration.

GR: coming soon!

https://github.com/user-attachments/assets/c85f0d53-7d4c-458e-8102-b37cc7d0abb7

## Install

```julia
] add https://github.com/JuliaAPlavin/Synchray.jl
```

## Minimal use (conical jet)

```julia
using Synchray
import Synchray as S
using RectiGrids

# Image-plane grid
xs = range(-1.0, 1.0; length=256)
ys = range(-1.0, 1.0; length=256)
cam = S.CameraZ(;
    xys=grid(S.SVector, xs, ys),
    nz=200, ν=230.0, t=0.0)

jet = S.EmissionRegion(
    geometry = S.Geometries.Conical(;
        axis = S.SVector(0.0, 0.0, 1.0),
        φj = 0.03,
        z = 0.1..10.0,
    ),
    velocity = S.VelocitySpec(S.Directions.Axial(), S.beta, S.Profiles.Constant(0.99)),
    emission = S.SynchrotronEmission(
        ne = S.Profiles.Axial(S.PowerLaw(-2; val0=1.0, s0=1.0)),
        B = S.BFieldSpec(
            S.Profiles.Axial(S.PowerLaw(-1; val0=1.0, s0=1.0)),
            S.Directions.Scalar(),
            b -> S.FullyTangled(b)
        ),
        electrons = S.IsotropicPowerLawElectrons(; p=2.5),
    ),
)

Iν = S.render(cam, S.prepare_for_computations(jet), S.Intensity())
# Iν is a 256×256 array of per-pixel specific intensities
```

`Unitful` inputs/outputs are supported: construct with unitful values then call `ustrip(region)` and `ustrip(cam)`, and use `S.withunits(S.render, cam, jet, S.Intensity())` to get unitful output intensity.

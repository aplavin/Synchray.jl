module Synchray

using Reexport
@reexport using StaticArrays
@reexport using IntervalSets
@reexport using Swizzling
@reexport using LinearAlgebra
using DispatchDoctor: @stable, @unstable

@stable default_mode="warn" begin
include("minkowski.jl")
include("mediums.jl")
include("objects.jl")
include("camera.jl")
include("transfer.jl")
end

export AbstractMedium,
	z_interval, four_velocity,
	emissivity, absorption,
	emissivity_invariant, absorption_invariant,
	UniformSphere,
		UniformSlab,
	OrthoCamera,
	integrate_ray,
	render,
	FourVector, photon_k, minkowski_dot, measured_frequency, doppler_factor

end

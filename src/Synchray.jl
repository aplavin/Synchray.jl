module Synchray

using Reexport
@reexport using StaticArrays
@reexport using IntervalSets
@reexport using Swizzling
@reexport using LinearAlgebra
@reexport using Unitful
@reexport using DataPipes
@reexport using Accessors
import AccessorsExtra: construct
using ForwardDiff
using Statistics: mean
using DispatchDoctor: @stable, @unstable


@stable default_mode="disable" begin
include("power.jl")
include("minkowski.jl")
include("mediums.jl")
include("camera.jl")
include("objects/slabs.jl")
include("objects/spheres.jl")
include("objects/bk_jet.jl")
include("transfer.jl")
end

# we don't export anything from this module on purpose
# the user is recommended to do `import Synchray as S`
# and then access things as `S.SomeType` or `S.some_function`

end

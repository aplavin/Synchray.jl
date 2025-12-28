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
include("camera.jl")
include("objects.jl")
include("transfer.jl")
end

# we don't export anything from this module on purpose
# the user is recommended to do `import Synchray as S`
# and then access things as `S.SomeType` or `S.some_function`

end

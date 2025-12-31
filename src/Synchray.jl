module Synchray

using Reexport
@reexport using StaticArrays
@reexport using IntervalSets
@reexport using Swizzling
@reexport using LinearAlgebra
@reexport using Unitful
using Unitful: AbstractQuantity
@reexport using DataPipes
@reexport using Accessors
import AccessorsExtra: construct
using ForwardDiff
import SpecialFunctions
using Statistics: mean
using StructArrays
using DispatchDoctor: @stable, @unstable


@stable default_mode="disable" begin
include("power.jl")
include("minkowski.jl")
include("mediums.jl")
include("camera.jl")
include("objects/slabs.jl")
include("objects/spheres.jl")
include("objects/bk_jet.jl")
include("patterns.jl")
include("transfer.jl")
include("units.jl")


to_float_type(::Type{TF}, obj) where {TF<:AbstractFloat} =
    @modify(x -> to_float_type(TF, x), obj[∗ₚ])
to_float_type(::Type{TF}, x::AbstractArray) where {TF<:AbstractFloat} =
    # XXX: can be more efficient for RectiGrids
    map(y -> to_float_type(TF, y), x)
to_float_type(::Type{TF}, x::AbstractFloat) where {TF<:AbstractFloat} = TF(x)

end

# we don't export anything from this module on purpose
# the user is recommended to do `import Synchray as S`
# and then access things as `S.SomeType` or `S.some_function`

end

module STOPIT

using StaticArrays
using ArgCheck
using QuadGK

include("units.jl")
include("particle.jl")
include("absorbers/absorbers.jl")
include("desorb/desorb.jl")

end # module STOPIT

module STOPIT

using StaticArrays
using ArgCheck
using QuadGK

include("units.jl")
include("particle.jl")
include("desorb/desorb.jl")

end # module STOPIT

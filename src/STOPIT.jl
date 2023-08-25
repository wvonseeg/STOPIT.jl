module STOPIT

using StaticArrays
using ArgCheck

include("units.jl")
include("errors.jl")
include("particle.jl")
include("absorbers/absorbers.jl")
include("desorb/desorb.jl")

end # module STOPIT

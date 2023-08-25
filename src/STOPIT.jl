module STOPIT

using StaticArrays
using ArgCheck

include("utils/units.jl")
include("utils/errors.jl")
include("particle.jl")
include("absorbers/absorbers.jl")
include("desorb/desorb.jl")

end # module STOPIT

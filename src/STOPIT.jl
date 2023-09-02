module STOPIT

using StaticArrays
using ArgCheck

include("utils/units.jl")
include("utils/errors.jl")
include("utils/masses.jl")
include("Particles/particles.jl")
include("Absorbers/absorbers.jl")
include("Interactive/interactive.jl")

using .Particles, .Absorbers

include("desorb.jl")

function __init__()
    Unitful.register(STOPIT)
end

end # module STOPIT
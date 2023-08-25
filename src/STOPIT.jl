module STOPIT

using StaticArrays
using ArgCheck
using QuadGK

include("units.jl")
include("errors.jl")
include("particle.jl")
include("absorbers/absorbers.jl")
include("desorb/desorb.jl")

# Added for testing for now
testsandwich = Sandwich()
silayer = getstandardmedium("Si"; depth=0.05u"cm")
proton = Particle(1, 1, 1.0, 10u"MeV")
setlayer!(testsandwich, silayer)

desorb!(proton, testsandwich)

end # module STOPIT

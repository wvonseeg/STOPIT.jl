module Particles
using NaturallyUnitful, ArgCheck
export Particle, info, Qvalue
struct Particle
    Z::UInt8
    A::UInt16
    mass::typeof(1.0u"u")
    energy::typeof(1.0u"MeV")
end

# Constructor to create particle by looking up mass from table
function Particle(Z::Integer, A::Integer, energy::typeof(1.0u"MeV"))
    mass = getmass(Z, A)
    Particle(Z, A, mass, energy)
end

function show(part::Particle)
    println("Z = ", part.Z)
    println("A = ", part.A)
    println("Mass = ", part.mass)
    println("E = ", part.energy)
end

include("reactions.jl")
end # module Particles
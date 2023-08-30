module Particles
using NaturallyUnitful, ArgCheck
import ..STOPIT: getmass
export Particle, info, Qvalue
struct Particle
    A::UInt16
    mass::typeof(1.0u"u")
    Z::UInt8
    energy::typeof(1.0u"MeV")
end

# Constructor to create particle by looking up mass from table
function Particle(A::Integer, Z::Integer, energy::typeof(1.0u"MeV"))
    mass = getmass(Z, A)
    Particle(Z, A, mass[1, 1], energy)
end

function show(part::Particle)
    println("Z = ", part.Z)
    println("A = ", part.A)
    println("Mass = ", part.mass)
    println("E = ", part.energy)
end

include("reactions.jl")
end # module Particles
module Particles
using NaturallyUnitful, ArgCheck
import ..STOPIT: getmass, atomicsymbols
export Particle, info, Qvalue, setenergy

struct Particle
    A::UInt16
    mass::typeof(1.0u"u")
    Z::UInt8
    energy::typeof(1.0u"MeV")
    symbol::String
end

# Constructor to create particle by looking up mass from table
function Particle(A::Integer, Z::Integer, energy::typeof(1.0u"MeV"))
    mass = getmass(Z, A)
    Particle(A, mass[1, 1], Z, energy, atomicsymbols[Z+1])
end

function Base.show(io::IO, ::MIME"text/plain", x::Particle)
    println("$(x.A)$(x.symbol) ($(x.mass)) @ $(x.energy)")
end

function setenergy(part::Particle, energy::typeof(1.0u"MeV"))
    return Particle(part.A, part.mass, part.Z, energy, part.symbol)
end

include("reactions.jl")
end # module Particles
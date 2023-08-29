export Particle, info
struct Particle
    Z::UInt8
    A::UInt16
    mass::Float64
    energy::typeof(1.0u"MeV")
end

# Constructor to create particle by looking up mass from table
function Particle(Z::Integer, A::Integer, energy::typeof(1.0u"MeV"))
    mass = getmass(Z, A)
    Particle(Z, A, mass, energy)
end

function info(part::Particle)
    print("Z = ", part.Z)
    print("A = ", part.A)
    print("E = ", part.energy)
end
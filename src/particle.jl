export Particle, info
struct Particle
    Z::UInt8
    Anumber::UInt16
    A::Float64
    energy::typeof(1.0u"MeV")
end

# Constructor to create particle by looking up mass from table

function info(part::Particle)
    print("Z = ", part.Z)
    print("A = ", part.A)
    print("E = ", part.energy)
end
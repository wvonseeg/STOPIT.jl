struct Particle
    Z::Integer
    Anumber::Integer
    A::Float64
    energy::typeof(1u"MeV")
end

# Constructor to create particle by looking up mass from table

function info(part::Particle)
    print("Z = ", part.Z)
    print("A = ", part.A)
    print("E = ", part.energy)
end
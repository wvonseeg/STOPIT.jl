module Particles
using NaturallyUnitful, ArgCheck
import ..STOPIT: getmass
export Particle, info, Qvalue

const atomicsymbols = ("n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
    "Lv", "Ts", "Og")

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
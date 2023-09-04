using NaturallyUnitful

import Unitful:
    nm, μm, mm, cm, m,
    μm, mg, g, kg, u,
    Ra, °F, °C, K,
    rad, °,
    ns, μs, ms, s, minute, hr, d,
    Hz, kHz, MHz,
    J, eV, keV, MeV, GeV,
    nA, μA, mA, A, mol,
    μW, mW, W, kW, C,
    mTorr, Torr, atm, bar, mbar, Pa, kPa

import Unitful:
    c0, mn, mp, q, me

import Unitful: 𝐋, 𝐌, 𝐓, 𝐈, 𝐉, 𝐍, 𝚯

import Unitful:
    Length, Area, Volume,
    Velocity, Luminosity, Density,
    Time, Frequency, Mass,
    Current, Charge, Power,
    Temperature, AbsoluteScaleTemperature, RelativeScaleTemperature

import Unitful: LengthUnits, AreaUnits, MassUnits, TemperatureUnits

const global atomicsymbols = ("n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs",
    "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
    "Lv", "Ts", "Og")

@derived_dimension Thickness 𝐌 / 𝐋^2 true
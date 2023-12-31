using Unitful

export masstoenergy

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
    c0, mn, mp, q

import Unitful: 𝐋, 𝐌, 𝐓, 𝐈, 𝐉, 𝐍, 𝚯

import Unitful:
    Length, Area, Volume,
    Velocity, Luminosity, Density,
    Time, Frequency, Mass,
    Current, Charge, Power,
    Temperature, AbsoluteScaleTemperature, RelativeScaleTemperature

import Unitful: LengthUnits, AreaUnits, MassUnits, TemperatureUnits

@derived_dimension Thickness 𝐌 / 𝐋^2 true

@inline function masstoenergy(mass::typeof(1.0u"u"))
    # 931.49410242(28) MeV
    return 931.49410242u"MeV/u" * mass
end
const global STANDARDMEDIA = ("C02", "Si", "Graphite", "C4H10", "CF4", "CH2", "CD2", "He", "H2", "Mylar", "P10")

"""
	getstandardmedium(name::String; pressure::Float64 = 0.0, depth::Float64 = 0.0)

Kwargs:
* `pressure` - Only need to input for gaseous targets. Value will be interpreted as Torr.
* `depth` - Necessary for both solid and gaseous targets. How long target is. Value will be interpreted as cm.

Available standard media are as follows:
* `C02` - Carbon dioxide...duh
* `Si` - Silicon, mostly so you don't have to remember the density
* `Graphite` - Just pure carbon, again for the density
* `C4H10` - Methane
* `CF4` - Not sure honestly
* `CH2` - Good targets
* `CD2` - Good, slightly heavier targets
* `He` - Gaseous helium, probably not entirely necessary
* `H2` - Gaseous hydrogen
* `Mylar` - Normally used for windows
* `P10` - Mix of 90% Argon and 10% Methane, useful in ionization chambers
"""
function getstandardmedium(name::String; pressure::Float64=0.0, depth::Float64=0.0)
    @argcheck name in STANDARDMEDIA
    depthwunits = depth * u"cm"
    pressurewunits = pressure * u"Torr"
    if name == "C02"
        return GasAbsorber((12, 16), (6, 8), (1, 2), pressurewunits, depthwunits)
    elseif name == "Si"
        dens = 2.329u"g/cm^3"
        return SolidAbsorber((28,), (14,), (1,), uconvert(u"mg/cm^2", dens * depthwunits), dens)
    elseif name == "Graphite"
        dens = 2.267u"g/cm^3"
        return SolidAbsorber((12,), (6,), (1,), uconvert(u"mg/cm^2", dens * depthwunits), dens)
    elseif name == "C4H10"
        return GasAbsorber((12, 1), (6, 1), (4, 10), pressurewunits, depthwunits)
    elseif name == "CF4"
        return GasAbsorber((12, 19), (6, 9), (1, 4), pressurewunits, depthwunits)
    elseif name == "CH2"
        dens = 0.94u"g/cm^3"
        return SolidAbsorber((12, 1), (6, 1), (1, 2), uconvert(u"mg/cm^2", dens * depthwunits), dens)
    elseif name == "CD2"
        dens = 0.94u"g/cm^3"
        return SolidAbsorber((12, 2), (6, 1), (1, 2), uconvert(u"mg/cm^2", dens * depthwunits), dens)
    elseif name == "He"
        return GasAbsorber((4,), (2,), (1,), pressurewunits, depthwunits)
    elseif name == "H2"
        return GasAbsorber((1,), (1,), (2,), pressurewunits, depthwunits)
    elseif name == "Mylar"
        dens = 1.38u"g/cm^3"
        return SolidAbsorber((12, 1, 16), (6, 1, 8), (10, 8, 4), uconvert(u"mg/cm^2", dens * depthwunits), dens)
    elseif name == "P10"
        return GasAbsorber((12, 1, 40), (6, 1, 18), (2, 8, 90), pressurewunits, depthwunits)
    end
end
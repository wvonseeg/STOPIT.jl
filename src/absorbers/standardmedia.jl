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
* `C4H10` - Butane, useful ionization gas
* `CF4` - Tetrafluoromethane aka carbon tetrafluoride, useful ionization gas
* `CH2` - Good proton target
* `CD2` - Good deuterium target
* `He` - Gaseous helium, probably not entirely necessary
* `H2` - Gaseous hydrogen
* `Mylar` - Normally used for windows
* `P10` - Mix of 90% Argon and 10% Methane, useful in ionization chambers
* `CH4` - Methane, pretty good ionization gas
"""
function getstandardmedium(name::String; pressure::Unitful.Pressure=0.0u"Torr", depth::Unitful.Length=0.0u"cm")
    @argcheck name in STANDARDMEDIA DomainError(name, "Not one of the available standard media")
    depth = uconvert(u"cm", depth)
    pressure = uconvert(u"Torr", pressure)
    if name == "C02"
        return GasAbsorber([12, 16], [6, 8], [1, 2], pressure, depth)
    elseif name == "Si"
        dens = 2.329u"g/cm^3"
        return SolidAbsorber([28,], [14,], [1,], dens * depth, dens)
    elseif name == "Graphite"
        dens = 2.267u"g/cm^3"
        return SolidAbsorber([12,], [6,], [1,], dens * depth, dens)
    elseif name == "C4H10"
        return GasAbsorber([12, 1], [6, 1], [4, 10], pressure, depth)
    elseif name == "CF4"
        return GasAbsorber([12, 19], [6, 9], [1, 4], pressure, depth)
    elseif name == "CH2"
        dens = 0.94u"g/cm^3"
        return SolidAbsorber([12, 1], [6, 1], [1, 2], dens * depth, dens)
    elseif name == "CD2"
        dens = 0.94u"g/cm^3"
        return SolidAbsorber([12, 2], [6, 1], [1, 2], dens * depth, dens)
    elseif name == "He"
        return GasAbsorber([4,], [2,], [1,], pressure, depth)
    elseif name == "H2"
        return GasAbsorber([1,], [1,], [2,], pressure, depth)
    elseif name == "Mylar"
        dens = 1.38u"g/cm^3"
        return SolidAbsorber([12, 1, 16], [6, 1, 8], [10, 8, 4], dens * depth, dens)
    elseif name == "P10"
        return GasAbsorber([12, 1, 40], [6, 1, 18], [2, 8, 90], pressure, depth)
    elseif name == "CH4"
        return GasAbsorber([12, 1], [6, 1], [1, 4], pressure, depth)
    end
end
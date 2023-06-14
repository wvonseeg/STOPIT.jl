const global STANDARDMEDIA = ("C02", "Si", "Graphite", "C4H10", "CF4", "CH2", "CD2", "He", "H2")

"""
	getstandardmedium(name::String; pressure::Float64 = 0.0, depth::Float64 = 0.0)


"""
function getstandardmedium(name::String; pressure::Float64=0.0, depth::Float64=0.0)
    @argcheck name in STANDARDMEDIA
    if name == "C02"
        return GasAbsorber((12, 16), (6, 8), (1, 2), pressure, depth)
    elseif name == "Si"
        return AbstractAbsorber((28,), (14,), (1,), 0.2329 * depth, 2.329)
    elseif name == "Graphite"
        return AbstractAbsorber((12,), (6,), (1,), 0.2267 * depth, 2.267)
    elseif name == "C4H10"
        return GasAbsorber((12, 1), (6, 1), (4, 10), pressure, depth)
    elseif name == "CF4"
        return GasAbsorber((12, 19), (6, 9), (1, 4), pressure, depth)
    elseif name == "CH2"
        return AbstractAbsorber((12, 1), (6, 1), (1, 2), 0.094 * depth, 0.94)
    elseif name == "CD2"
        return AbstractAbsorber((12, 2), (6, 1), (1, 2), 0.094 * depth, 0.94)
    elseif name == "He"
        return GasLayer((4,), (2,), (1,), pressure, depth)
    elseif name == "H2"
        return GasLayer((1,), (1,), (2,), pressure, depth)
    end
end
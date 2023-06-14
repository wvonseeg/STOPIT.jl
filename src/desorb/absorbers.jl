export AbsorberLayer, GasAbsorber, Sandwich, setlayer!

struct AbsorberLayer
    A::Tuple{Float64}
    Z::Tuple{Integer}
    num::Tuple{Integer}
    thickness::Float64
    partialthickness::Tuple{Float64}
    density::Float64
    partialdensity::Tuple{Float64}
end

function AbsorberLayer(A::Vector{Real}, Z::Vector{Integer}, num::Vector{Real}, thickness::Float64, density::Float64)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)

    len = length(A)
    ANout = zeros(len)
    Tout = zeros(len)
    Dout = zeros(len)

    for i = 1:len
        ANout[i] = A[i] * num[i]
        Tout[i] = thickness * num[i]
        Dout[i] = density * num[i]
    end

    ANout ./= sum(ANout)
    AbsorberLayer(tuple(A), tuple(Z), tuple(ANout), thickness, tuple(Tout), density, tuple(Dout))
end

struct GasAbsorber <: AbsorberLayer
    A::Tuple{Float64}
    Z::Tuple{Integer}
    num::Tuple{Integer}
    thickness::Float64
    partialthickness::Tuple{Float64}
    density::Float64
    partialdensity::Tuple{Float64}
    concentrations::Tuple{Float64}
    pressure::Float64
    depth::Float64
end

function GasAbsorber(A::Vector{Real}, Z::Vector{Integer}, num::Vector{Real}, concentrations::Vector{Real}, pressure::Float64, depth::Float64)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)
    @argcheck length(concentrations) == length(num)

    len = length(A)
    thick = zeros(len)
    dens = zeros(len)

    P = PR / 760.0
    X = XL / 22.4

    # AWW = 0.0
    # AW = 0.0
    thickness = 0.0
    density = 0.0
    for i = 1:INW
        # AW += A[i] * AN[i]
        # AWW += A[i] * AN[i] * CN[i]
        thick[i] = P * X * A[i] * AN[i] * CN[i]
        thickness += thick[i]
        dens[i] = thick[i] / XL
        density += dens[i]
    end
    GasAbsorber(tuple(A), tuple(Z), tuple(num), thickness, tuple(thick), density, tuple(dens), tuple(concentrations), pressure, depth)
end

function GasAbsorber(A::Vector{Real}, Z::Vector{Integer}, num::Vector{Real}, pressure::Float64, depth::Float64)
    concentrations = ones(length(A))
    GasAbsorber(A, Z, num, concentrations, pressure, depth)
end

mutable struct Sandwich
    layers::Vector{AbsorberLayer}
end

function setlayer!(sandwich::Sandwich, A::Vector{Real}, Z::Vector{Integer}, AN::Vector{Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, AbsorberLayer(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, A::Vector{Real}, Z::Vector{Integer}, AN::Vector{Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    @argcheck 0 < index <= length(sandwich.layers) + 1
    insert!(sandwich.layers, index, AbsorberLayer(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, A::Vector{Real}, Z::Vector{Integer}, AN::Vector{Integer}, CN::Vector{Real}, pressure::Float64, depth::Float64)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, GasAbsorber(A, Z, AN, CN, pressure, depth))
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, A::Vector{Real}, Z::Vector{Integer}, AN::Vector{Integer}, CN::Vector{Real}, pressure::Float64, depth::Float64)
    @argcheck length(sandwich.layers) < 19
    @argcheck 0 < index <= length(sandwich.layers) + 1
    insert!(sandwich.layers, index, GasAbsorber(A, Z, AN, CN, pressure, depth))
    return sandwich
end

function removelayer!(sandwich::Sandwich)
    @argcheck length(sandwich.layers) > 0
    pop!(sandwich.layers)
    return sandwich
end

function removelayer!(sandwich::Sandwich, index::Integer)
    @argcheck length(sandwich.layers) > 0
    @argcheck 0 < index <= length(sandwich.layers)
    popat!(sandwich.layers, index)
    return sandwich
end
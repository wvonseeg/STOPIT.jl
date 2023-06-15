export AbstractAbsorber, SolidAbsorber, GasAbsorber, Sandwich, setlayer!, removelayer!

"""
Abstract type for absorber layers.
"""
abstract type AbstractAbsorber end

"""
Subtype of AbstractAbsorber for solid absorber layers.

The units of fields are as follows:
* Thicnkess => mg cm^-2
* Density   => g cm^-3
"""
struct SolidAbsorber <: AbstractAbsorber
    A::Tuple{Float64}
    Z::Tuple{Integer}
    num::Tuple{Integer}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Tuple{typeof(1.0u"mg/cm^2")}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Tuple{typeof(1.0u"g/cm^3")}
end

function SolidAbsorber(A::Vector{Float64}, Z::Vector{Integer}, num::Vector{Float64}, thickness::Float64, density::Float64)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)

    len = length(A)
    ANout = zeros(len)
    Tout = zeros(typeof(1.0u"mg/cm^2"), len)
    Dout = zeros(typeof(1.0u"g/cm^3"), len)

    for i = 1:len
        @inbounds ANout[i] = A[i] * num[i]
        @inbounds Tout[i] = thickness * num[i] * 1u"mg/cm^2"
        @inbounds Dout[i] = density * num[i] * 1u"g/cm^2"
    end

    ANout ./= sum(ANout)
    SolidAbsorber(tuple(A), tuple(Z), tuple(ANout), thickness * 1u"mg/cm^2", tuple(Tout), density * 1u"g/cm^3", tuple(Dout))
end

"""
Subtype of AbstractAbsorber for gaseous absorber layers.

The units of fields are as follows:
* Thickness => mg cm^-2
* Density   => g cm^-3
* Pressure  => Torr
* Depth     => cm
"""
struct GasAbsorber <: AbstractAbsorber
    A::Tuple{Float64}
    Z::Tuple{Integer}
    num::Tuple{Integer}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Tuple{typeof(1.0u"mg/cm^2")}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Tuple{typeof(1.0u"g/cm^3")}
    concentrations::Tuple{Float64}
    pressure::typeof(1.0u"Torr")
    depth::typeof(1.0u"cm")
end

function GasAbsorber(A::Vector{Float64}, Z::Vector{Integer}, num::Vector{Integer}, concentrations::Vector{Float64}, pressure::Float64, depth::Float64)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)
    @argcheck length(concentrations) == length(num)

    len = length(A)
    thick = zeros(typeof(1.0u"mg/cm^2"), len)
    dens = zeros(typeof(1.0u"g/cm^3"), len)

    P = pressure / 760.0
    X = depth / 22.4

    # AWW = 0.0
    # AW = 0.0
    thickness = 0.0u"mg/cm^2"
    density = 0.0u"g/cm^3"
    for i = 1:len
        # AW += A[i] * AN[i]
        # AWW += A[i] * AN[i] * CN[i]
        @inbounds thick[i] = P * X * A[i] * num[i] * concentrations[i] * 1u"mg/cm^2"
        @inbounds thickness += thick[i]
        @inbounds dens[i] = thick[i] / depth * 1u"g/cm^3"
        @inbounds density += dens[i]
    end
    GasAbsorber(tuple(A), tuple(Z), tuple(num), thickness, tuple(thick), density, tuple(dens), tuple(concentrations), pressure * 1u"Torr", depth * 1u"cm")
end

function GasAbsorber(A::Vector{Float64}, Z::Vector{Integer}, num::Vector{Integer}, pressure::Float64, depth::Float64)
    concentrations = ones(length(A))
    GasAbsorber(A, Z, num, concentrations, pressure, depth)
end

mutable struct Sandwich
    layers::Vector{AbstractAbsorber}
end

function setlayer!(sandwich::Sandwich, A::Vector{Float64}, Z::Vector{Integer}, AN::Vector{Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, SolidAbsorber(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, A::Vector{Float64}, Z::Vector{Integer}, AN::Vector{Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    @argcheck 0 < index <= length(sandwich.layers) + 1
    insert!(sandwich.layers, index, SolidAbsorber(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, A::Vector{Float64}, Z::Vector{Integer}, AN::Vector{Integer}, CN::Vector{Float64}, pressure::Float64, depth::Float64)
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
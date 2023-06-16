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
    Z::Tuple{UInt8}
    num::Tuple{UInt8}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Tuple{typeof(1.0u"mg/cm^2")}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Tuple{typeof(1.0u"g/cm^3")}
end

#################################
### SolidAbsorber Constructors###
#################################

function SolidAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{<:Integer}, thickness::Thickness, density::Unitful.Density)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)

    thickness = uconvert(u"mg/cm^2", thickness)
    density = uconvert(u"g/cm^3", density)

    len = length(A)
    ANout = zeros(len)
    Tout = zeros(typeof(1.0u"mg/cm^2"), len)
    Dout = zeros(typeof(1.0u"g/cm^3"), len)

    @inbounds for i = 1:len
        ANout[i] = A[i] * num[i]
        Tout[i] = thickness * num[i]
        Dout[i] = density * num[i]
    end

    ANout ./= sum(ANout)
    SolidAbsorber(tuple(A...), tuple(convert(Vector{UInt8}, Z)...),
        tuple(convert(Vector{UInt8}, ANout)...), thickness,
        tuple(Tout...), density, tuple(Dout...))
end

function SolidAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, thickness::Thickness, density::Unitful.Density)
    SolidAbsorber(convert(Vector{Float64}, A), Z, num, thickness, density)
end

#function SolidAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{Float64}, thickness::Float64, density::Float64)
#    SolidAbsorber(A, Z, num, thickness * 1.0u"mg/cm^2", density * 1.0u"g/cm^3")
#end

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
    Z::Tuple{UInt8}
    num::Tuple{UInt8}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Tuple{typeof(1.0u"mg/cm^2")}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Tuple{typeof(1.0u"g/cm^3")}
    concentrations::Tuple{Float64}
    pressure::typeof(1.0u"Torr")
    depth::typeof(1.0u"cm")
end

#################################
#### GasAbsorber Constructors####
#################################

function GasAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{<:Integer}, concentrations::Vector{Float64}, pressure::Unitful.Pressure, depth::Unitful.Length)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    @argcheck length(A) <= 4
    @argcheck length(Z) == length(A)
    @argcheck length(num) == length(Z)
    @argcheck length(concentrations) == length(num)

    pressure = uconvert(u"Torr", pressure)
    depth = uconvert(u"cm", depth)


    len = length(A)
    thick = zeros(typeof(1.0u"mg/cm^2"), len)
    dens = zeros(typeof(1.0u"g/cm^3"), len)

    P = pressure / 760.0u"Torr"
    X = depth / 22.4u"cm"

    # AWW = 0.0
    # AW = 0.0
    thickness = 0.0u"mg/cm^2"
    density = 0.0u"g/cm^3"
    @inbounds for i = 1:len
        # AW += A[i] * AN[i]
        # AWW += A[i] * AN[i] * CN[i]
        thick[i] = P * X * A[i] * num[i] * concentrations[i] * 1u"mg/cm^2"
        thickness += thick[i]
        dens[i] = thick[i] / depth
        density += dens[i]
    end
    GasAbsorber(tuple(A...), tuple(convert(Vector{UInt8}, Z)...),
        tuple(convert(Vector{UInt8}, num)...), thickness, tuple(thick...),
        density, tuple(dens...), tuple(concentrations...), pressure, depth)
end

function GasAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{<:Integer}, concentrations::Vector{Float64}, pressure::Float64, depth::Float64)
    GasAbsorber(A, Z, num, concentrations, pressure * 1u"Torr", depth * 1u"cm")
end

function GasAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{<:Integer}, pressure::Unitful.Pressure, depth::Unitful.Length)
    concentrations = ones(length(A))
    GasAbsorber(A, Z, num, concentrations, pressure, depth)
end

function GasAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, pressure::Unitful.Pressure, depth::Unitful.Length)
    GasAbsorber(convert(Vector{Float64}, A), Z, num, pressure, depth)
end

function GasAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{<:Integer}, pressure::Float64, depth::Float64)
    GasAbsorber(A, Z, num, pressure * 1u"Torr", depth * 1u"cm")
end

function GasAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, pressure::Float64, depth::Float64)
    GasAbsorber(convert(Vector{Float64}, A), Z, num, pressure, depth)
end
mutable struct Sandwich
    layers::Vector{<:AbstractAbsorber}
end

function Sandwich()
    Sandwich(Vector{AbstractAbsorber}())
end

function setlayer!(sandwich::Sandwich, layer::AbstractAbsorber)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, layer)
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, layer::AbstractAbsorber)
    @argcheck length(sandwich.layers) < 19
    @argcheck 0 < index <= length(sandwich.layers) + 1
    insert!(sandwich.layers, index, layer)
    return sandwich
end

function setlayer!(sandwich::Sandwich, A::Vector{Float64}, Z::Vector{<:Integer}, AN::Vector{<:Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, SolidAbsorber(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, A::Vector{Float64}, Z::Vector{<:Integer}, AN::Vector{<:Integer}, thickness::Float64, density::Float64)
    @argcheck length(sandwich.layers) < 19
    @argcheck 0 < index <= length(sandwich.layers) + 1
    insert!(sandwich.layers, index, SolidAbsorber(A, Z, AN, thickness, density))
    return sandwich
end

function setlayer!(sandwich::Sandwich, A::Vector{Float64}, Z::Vector{<:Integer}, AN::Vector{<:Integer}, CN::Vector{Float64}, pressure::Float64, depth::Float64)
    @argcheck length(sandwich.layers) < 19
    push!(sandwich.layers, GasAbsorber(A, Z, AN, CN, pressure, depth))
    return sandwich
end

function setlayer!(sandwich::Sandwich, index::Integer, A::Vector{Float64}, Z::Vector{<:Integer}, AN::Vector{<:Integer}, CN::Vector{Float64}, pressure::Float64, depth::Float64)
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
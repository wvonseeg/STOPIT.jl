module Absorbers
export AbstractAbsorber,
    SolidAbsorber,
    GasAbsorber,
    Sandwich,
    setlayer!,
    removelayer!,
    purgelayers!,
    getstandardmedium

"""
Abstract type for absorber layers.
"""
abstract type AbstractAbsorber end

mutable struct Sandwich
    layers::Vector{<:AbstractAbsorber}
end

function Sandwich()
    Sandwich(Vector{AbstractAbsorber}())
end

include("solids.jl")
include("gases.jl")
include("standardmedia.jl")

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

function purgelayers!(sandwich::Sandwich)
    empty!(sandwich.layers)
    return sandwich
end
end # module Absorbers
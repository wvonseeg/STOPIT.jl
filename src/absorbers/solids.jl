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
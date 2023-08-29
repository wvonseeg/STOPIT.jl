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
    len = length(A)
    @argcheck len <= 4
    @argcheck length(Z) == len DimensionMismatch
    @argcheck length(num) == len DimensionMismatch

    thickness = uconvert(u"mg/cm^2", thickness)
    density = uconvert(u"g/cm^3", density)

    ANout = zeros(len)
    Tout = zeros(typeof(1.0u"mg/cm^2"), len)
    Dout = zeros(typeof(1.0u"g/cm^3"), len)

    @inbounds for (i, N) in enumerate(num)
        ANout[i] = A[i] * N
        Tout[i] = thickness * N
        Dout[i] = density * N
    end

    ANout ./= sum(ANout)
    SolidAbsorber(tuple(A...), tuple(convert(Vector{UInt8}, Z)...),
        tuple(convert(Vector{UInt8}, ANout)...), thickness,
        tuple(Tout...), density, tuple(Dout...))
end

function SolidAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, thickness::Thickness, density::Unitful.Density)
    massmatrix = getmass(A, Z)
    SolidAbsorber(massmatrix[:, 1], Z, num, thickness, density)
end

#function SolidAbsorber(A::Vector{Float64}, Z::Vector{<:Integer}, num::Vector{Float64}, thickness::Float64, density::Float64)
#    SolidAbsorber(A, Z, num, thickness * 1.0u"mg/cm^2", density * 1.0u"g/cm^3")
#end
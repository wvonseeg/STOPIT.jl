"""
Subtype of AbstractAbsorber for solid absorber layers.

The units of fields are as follows:
* Thicnkess => mg cm^-2
* Density   => g cm^-3
"""
struct SolidAbsorber <: AbstractAbsorber
    A::Tuple{UInt16}
    mass::Tuple{typeof(1.0u"u")}
    Z::Tuple{UInt8}
    num::Tuple{UInt8}
    symbols::Tuple{String}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Tuple{typeof(1.0u"mg/cm^2")}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Tuple{typeof(1.0u"g/cm^3")}
end

#################################
### SolidAbsorber Constructors###
#################################

function SolidAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, thickness::typeof(1.0u"mg/cm^2"), density::Unitful.Density)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    len = length(A)
    @argcheck len <= 4
    @argcheck length(Z) == len DimensionMismatch
    @argcheck length(num) == len DimensionMismatch

    #thickness = uconvert(u"mg/cm^2", thickness)
    density = uconvert(u"g/cm^3", density)

    ANout = zeros(len)
    Tout = zeros(typeof(1.0u"mg/cm^2"), len)
    Dout = zeros(typeof(1.0u"g/cm^3"), len)

    syms = Vector{String}(undef, len)

    @inbounds for (i, N) in enumerate(num)
        ANout[i] = A[i] * N
        Tout[i] = thickness * N
        Dout[i] = density * N
        syms[i] = atomicsymbols[Z[i]+1]
    end

    ANout ./= sum(ANout)
    massmatrix = getmass(A, Z)
    SolidAbsorber(tuple(A...), tuple(massmatrix[:, 1]...),
        tuple(convert(Vector{UInt8}, Z)...),
        tuple(convert(Vector{UInt8}, ANout)...), tuple(syms...),
        thickness, tuple(Tout...), density, tuple(Dout...))
end

function SolidAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, depth::Unitful.Length, density::Unitful.Density)
    SolidAbsorber(A, Z, num, uconvert(u"mg/cm^2", density * depth), density)
end

function SolidAbsorber(A::Vector{<:Real}, Z::Vector{<:Real}, num::Vector{<:Real}, thickness::Real, density::Real)
    SolidAbsorber(convert(Vector{Int}, A), convert(Vector{Int}, Z),
        convert(Vector{Int}, num), convert(Float64, thickness) * 1.0u"mg/cm^2", convert(Float64, density) * 1.0u"g/cm^3")
end

function Base.show(io::IO, ::MIME"text/plain", x::SolidAbsorber)
    println("Solid Absorber")
    println("Density: $(x.density)\tThickness: $(x.thickness)")
    println("Molecular Makeup")
    for (i, sym) in enumerate(x.symbols)
        println("$(x.num[i])\t$(x.A[i])$sym")
    end
end
"""
Subtype of AbstractAbsorber for solid absorber layers.

The units of fields are as follows:
* Thicnkess => mg cm^-2
* Density   => g cm^-3
"""
struct SolidAbsorber <: AbstractAbsorber
    A::Array{UInt16,1}
    mass::Array{typeof(1.0u"u"),1}
    Z::Array{UInt8,1}
    num::Array{UInt8,1}
    symbols::Array{String,1}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Array{typeof(1.0u"mg/cm^2"),1}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Array{typeof(1.0u"g/cm^3"),1}
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

    Tout = zeros(typeof(1.0u"mg/cm^2"), len)
    Dout = zeros(typeof(1.0u"g/cm^3"), len)

    syms = Vector{String}(undef, len)

    @inbounds for (i, N) in enumerate(num)
        Tout[i] = thickness * N
        Dout[i] = density * N
        syms[i] = atomicsymbols[Z[i]+1]
    end

    Tout ./= sum(num)
    Dout ./= sum(num)

    massmatrix = getmass(A, Z)
    SolidAbsorber(A, massmatrix[:, 1], Z, num, syms, thickness, Tout, density, Dout)
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

function setthickness(absorb::SolidAbsorber, thickness::typeof(1.0u"mg/cm^2"))
    pthick = Vector{typeof(1.0u"mg/cm^2")}(undef, length(absorb.partialthickness))
    for (i, N) in enumerate(absorb.num)
        pthick[i] = absorb.thickness * N
    end
    pthick ./= sum(num)
    return SolidAbsorber(absorb.A, absorb.mass, absorb.Z, absorb.num,
        absorb.symbols, thickness, pthick, absorb.density, absorb.partialdensity)
end
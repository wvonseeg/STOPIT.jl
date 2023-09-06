"""
Subtype of AbstractAbsorber for gaseous absorber layers.

The units of fields are as follows:
* Thickness => mg cm^-2
* Density   => g cm^-3
* Pressure  => Torr
* Depth     => cm
"""
struct GasAbsorber <: AbstractAbsorber
    A::Array{UInt16,1}
    mass::Array{typeof(1.0u"u"),1}
    Z::Array{UInt8,1}
    num::Array{UInt8,1}
    symbols::Array{String,1}
    thickness::typeof(1.0u"mg/cm^2")
    partialthickness::Array{typeof(1.0u"mg/cm^2"),1}
    density::typeof(1.0u"g/cm^3")
    partialdensity::Array{typeof(1.0u"g/cm^3"),1}
    concentrations::Array{Float64,1}
    pressure::typeof(1.0u"Torr")
    depth::typeof(1.0u"cm")
end

#################################
#### GasAbsorber Constructors####
#################################

function GasAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, concentrations::Vector{Float64}, pressure::Unitful.Pressure, depth::Unitful.Length)
    # Check arguments
    # For now only allow up to 4 elements in composite absorber
    # This number is taken as legacy from FROTRAN code and may be modified later
    len = length(A)
    @argcheck len <= 4
    @argcheck length(Z) == len DimensionMismatch
    @argcheck length(num) == len DimensionMismatch
    @argcheck length(concentrations) == len DimensionMismatch

    pressure = uconvert(u"Torr", pressure)
    depth = uconvert(u"cm", depth)

    thick = zeros(typeof(1.0u"mg/cm^2"), len)
    dens = zeros(typeof(1.0u"g/cm^3"), len)

    P = pressure / 760.0u"Torr"
    X = depth / 22.4u"cm"

    # AWW = 0.0
    # AW = 0.0
    thickness = 0.0u"mg/cm^2"
    density = 0.0u"g/cm^3"

    syms = Vector{String}(undef, len)
    @inbounds for i = 1:len
        # AW += A[i] * AN[i]
        # AWW += A[i] * AN[i] * CN[i]
        thick[i] = P * X * A[i] * num[i] * concentrations[i] * 1u"mg/cm^2"
        thickness += thick[i]
        dens[i] = thick[i] / depth
        density += dens[i]
        syms[i] = atomicsymbols[Z[i]+1]
    end

    thick ./= sum(num)
    dens ./= sum(num)

    massmatrix = getmass(A, Z)
    GasAbsorber(A, massmatrix[:, 1], Z, num, syms, thickness, thick,
        density, dens, concentrations, pressure, depth)
end

function GasAbsorber(A::Vector{<:Real}, Z::Vector{<:Real}, num::Vector{<:Real}, concentrations::Vector{Real}, pressure::Real, depth::Real)
    GasAbsorber(A, Z, num, concentrations, pressure * 1u"Torr", depth * 1u"cm")
end

function GasAbsorber(A::Vector{<:Integer}, Z::Vector{<:Integer}, num::Vector{<:Integer}, pressure::Unitful.Pressure, depth::Unitful.Length)
    concentrations = ones(length(A))
    GasAbsorber(A, Z, num, concentrations, pressure, depth)
end

function GasAbsorber(A::Vector{<:Real}, Z::Vector{<:Real}, num::Vector{<:Real}, pressure::Real, depth::Real)
    concentrations = ones(length(A))
    GasAbsorber(A, Z, num, concentrations, pressure, depth)
end

function Base.show(io::IO, ::MIME"text/plain", x::GasAbsorber)
    println("Gaseous Absorber")
    println("Pressure: $(x.pressure)\tDepth: $(x.depth)")
    for (i, sym) in enumerate(x.symbols)
        print(x.num[i])
        println("\t$(x.A[i])$sym")
    end
end
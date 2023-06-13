export CompositeAbsorber, SolidSandwich, GasSandwich, setabs, setabg

struct CompositeAbsorber
    A::Tuple{Float64}
    Z::Tuple{Integer}
    num::Tuple{Integer}
    thickness::Float64
    thick::Tuple{Float64}
    density::Float64
    dens::Tuple{Float64}
end

function CompositeAbsorber(A::Vector{Real}, Z::Vector{Integer}, num::Vector{Real}, thickness::Float64, density::Float64)
    ANout = A * num
    ANout /= sum(ANout)
    Tout = thickness * num
    Dout = density * num
    CompositeAbsorber(tuple(A), tuple(Z), tuple(num), thickness, tuple(Tout), density, tuple(Dout))
end

mutable struct SolidSandwich

end

mutable struct GasSandwich

end

function setabs(INW, A, Z, AN, TH, DN)
    # 	SUBROUTINE FOR SETTING UP COMPOSITE ABSORBER
    # 	DATA (PARTIAL DENSITIES AND THICKNESSES)

    #	DIMENSION A(4),Z(4),AN(4),T(4),D(4)
    AW = 0.0

    for i = 1:INW
        AW += A[i] * AN[i]
        AN[i] = A[i] * AN[i]
        T[i] = TH * AN[i]
        D[i] = DN * AN[i]
    end
    AN ./ AW
end

function setabg(INW, A, Z, AN, CN, TH, DN, PR, XL)
    # 	SUBROUTINE FOR SETTING UP COMPOSITE ABSORBER DATA
    # 	FOR GASEOUS LAYERS.

    #	DIMENSION A(4),Z(4),AN(4),CN(4),T(4),D(4)
    temp = []
    dens = []

    P = PR / 760.0
    X = XL / 22.4

    AWW = 0.0
    AW = 0.0

    for i = 1:INW
        AW += A[i] * AN[i]
        AWW += A[i] * AN[i] * CN[i]
        temp[i] = P * X * A[i] * AN[i] * CN[i]
        dens[i] = temp[i] / XL
    end
    return temp, dens
end
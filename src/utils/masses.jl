using Artifacts

export getmass

const AMEpath = joinpath(artifact"AME2020", "ame2020.txt")
const AZpattern = r"^(?:0+|\s+)\s+(?:\d+|-\d+)\s+\d+\s+(\d+)\s+(\d+)"
const masspattern = r"(\d+)\s(\d+\.\d+|\d+(?=#))"
const sigmapattern = r"(\d+\.\d+$|\d+(?=#$))"

function getmass(A::Integer, Z::Integer)
    return getmass([Z], [A])
end

function getmass(A::Vector{<:Integer}, Z::Vector{<:Integer})
    @argcheck length(Z) == length(A) DimensionMismatch

    order = sortperm(A)
    Zsearch = Z[order]
    Asearch = A[order]

    AZcounter = 1
    masslist = Array{Float64,2}(undef, length(Z), 2)
    open(AMEpath) do massfile
        # Read through each line looking for the correct Z & A combo
        for line in eachline(massfile)
            AZmatch = match(AZpattern, line)
            if AZmatch !== nothing
                Zmatch = parse(Int, AZmatch.captures[1])
                Amatch = parse(Int, AZmatch.captures[2])
                # If A & Z match, then extract mass and uncertainty
                if Amatch == Asearch[AZcounter] && Zmatch == Zsearch[AZcounter]
                    masslist[AZcounter, 1], masslist[AZcounter, 2] = _extractmasses(line)
                    if AZcounter == length(Asearch)
                        return masslist
                    end
                    AZcounter += 1
                elseif Amatch > Asearch[AZcounter]
                    # If we have searched past A, then the isotope is not in the file
                    throw(MassNotFoundException())
                end
            end
        end
    end
end

function _extractmasses(AMEline::String)
    massmatch = match(masspattern, AMEline)
    sigmamatch = match(sigmapattern, AMEline)
    if massmatch !== nothing || sigmamatch !== nothing
        mass = parse(Float64, massmatch.captures[1])
        mass += endswith(massmatch.captures[2], "#") ? 0.0 : parse(Float64, massmatch.captures[2]) * 1e-6
        sigma = endswith(sigmamatch.captures[1], "#") ? 0.0 : parse(Float64, sigmamatch.captures[1]) * 1e-6
        return mass, sigma
    end
end
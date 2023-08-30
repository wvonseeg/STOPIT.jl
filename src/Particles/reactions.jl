const RXNTYPES = ("g", "nn", "p", "a")

function Qvalue(target::Particle, projectile::Particle, rxntype::String)
    productA = target.A + projectile.A
    productZ = target.Z + projectile.Z
    productmass = 0.0
    if rxntype == "g"
        productmass = getmass(productA, productZ)[1, 1]
    elseif rxntype == "nn"
        productA -= 1
        productmass = getmass([productA,1], [productZ,0])[1, 1]
    elseif rxntype == "p"
        productA -= 1
        productZ -= 1
        productmass = getmass([productA,1], [productZ,1])[1, 1]
    elseif rxntype == "a"
        productA -= 4
        productZ -= 2
        productmass = getmass([productA,4], [productZ,2])[1, 1]
    end
    return masstoenergy(target.mass + projectile.mass - productmass)
end
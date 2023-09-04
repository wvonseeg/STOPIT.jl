export desorb!, desorb, dedx, ads, stoppingenergy, stoppingdepth

function desorb!(lostenergy::Vector{typeof(1.0u"MeV")}, stopee::Particle, sandwich::Sandwich;
    integrationsteps::Integer=1000, precision::Float64=1E-4, inplace=false)

    tmpstopee = stopee

    @inbounds for i = eachindex(sandwich.layers)
        try
            lostenergy[i] = ads(tmpstopee, sandwich.layers[i], integrationsteps, precision)
        catch err
            if isa(err, ParticleStoppedException)
                lostenergy[i] = tmpstopee.energy
                for j in i:length(sandwich.layers)
                    lostenergy[i] = 0.0u"MeV"
                end
                tmpstopee = setenergy(tmpstopee, tmpstopee.energy - lostenergy[i])
                break
            else
                rethrow
            end
        end
        tmpstopee = setenergy(tmpstopee, tmpstopee.energy - lostenergy[i])
    end

    if inplace
        stopee = tmpstopee
    end
    return lostenergy
end

function desorb(stopee::Particle, sandwich::Sandwich;
    integrationsteps::Integer=1000, precision::Float64=1e-4)
    lostenergy = Vector{typeof(1.0u"MeV")}(undef, length(sandwich.layers))
    desorb!(lostenergy, stopee, sandwich; integrationsteps=integrationsteps, precision=precision)
    return lostenergy
end

#=
ORIGINAL FORTRAN CODE TRANSLATION
function desorb(ianz, zp, ap, ep, loste)

    # **  CALCULATES ENERGY LOSS IN AN ABSORBER SANDWICH *******
    # ****  ENERGY DEPOSIT IN SECTIONS OF IONIZATION  **********
    # ***************  CHAMBER (DE1, DE2, AND DE3)  ************
    #
    #	save eltimo,izmax,izth

    #       the mass table is to be used only for iopt = 5,6
    #       use atomic masses to average for isotipic composition.
    #       taken from Formulas Facts and Constants, H. J. Fischbeck and
    #       K. H. Fischbeck. Springer - Verlag 1987 2nd ed, pages 164-183.


    #      IOPT = 1 - SUPPLY ENERGY OF PARTICLE ENTERING
    #                 THE ABSORBER ARRAY AND GET LOSS AND
    #                 RANGES
    #      IOPT = 2 - SUPPLY TARGET, PROECTILE AND EJECTILE
    #                 INFO. AND THEN GO THROUGH ABSORBER
    #                 SANDWICH
    #      IOPT = 3 - CALCULATE ENERGY DEPOSITS IN DETECTOR
    #                 DETECTOR DIMENSIONS ARE STANDARD AND
    #                 THE VARIABLE -'IDET' - CHOOSES BETWEEN
    #                 VARIETY OF AVAILABLE DETECTORS
    #      IOPT = 4 - FINDS MAXIMUM ENERGY THAT CAN BE STOPPED IN
    #                 IANZ ELEMENTS OF THE SANDWICH FOR GIVEN
    #                 ZP, AP.
    #                 WHEN CALCULATION IS FINISHED, THE PROGRAM READS
    #                 IN NEW VALUES OF ZP, AP AND RESTARTS. TO END
    #                 THE PROGRAM, GIVE ZP < 0.
    #                 IN ORDER TO HELP THE SPEED OF THE PROGRAM,
    #                 GIVE THE PARTICLE'S "Z" IN increasing  ORDER.
    #      IOPT = 5 - STORES ARRAYS OF Edet AS A FUNCTION OF INCIDENT    
    #                 ENERGY AND THE PARTICLE'S ID (Z,A)
    #                 ARRAY LABELED  eptable(Z-Zth,Einc,ipunch)
    #                 ipunch = 1 stopped,  = 2 punched through
    #                 Einc = E(incident)/detable
    #                 Zth  = lowest Z considered - 1
    #
    # ************************************************************
    #
    @argcheck 0 < iopt <= 6


    # ***********************************************************
    #
    #	if(iopt.eq.3) then
    #      open(unit=IO3, status='OLD', file='absorbv.out')
    #      rewind IO3
    #	endif
    #
    #	if(iopt.ge.5) then
    #		do itpun = 1,2
    #		do itblz = 1,50
    #		emintabz(itblz) = 0.
    #		do itble = 1,500
    #		eptable ( itblz,itble,itpun ) = 0.
    #		enddo
    #		enddo
    #		enddo
    #
    #		open(unit=io0, status='UNKNOWN', file='abs.tbl')
    #		rewind io0
    #
    #	else
    #	endif

    # ***********************************************************

    #      ianz  = number of elements in absorber "sandwich" - the 
    #              particle deposits all its energy in these layers.
    #      ianzi = index of last layer in which the energy of the 
    #              particle is not recorded - this unreecorded energy
    #              is used in the DST production in two modes:
    #          when making DST tape from data:
    #              Since the detector records only deposited energy
    #              the output tables are used to correct for this 
    #              deficinecy and for a given dharge and mass extrapolate
    #              the measured energy to the target exit energy
    #          when making DST tape from model calculations:
    #              The "lost" energy is a source of broadening of the 
    #              energy spectra due to straggling - this "smearing" 
    #              is estimated and superposed on the calculated spectra.
    #	ianzide = element # for DE calculation
    #
    #	ianzv(5)
    #	     = Index of layer where exiting particle velocit is 
    #	       calculated (only for option 3) max = 5
    # ***********************************************************


    #      AP = PROJECTILE MASS
    #      ZP = PROJECTILE CHARGE
    #      EP = PROJECTILE ENERGY

    # *************************************************************
    #
    #      zth= threshold Z for incident energy tble calc. (iopt=5,6)
    #      zmax = maximum z for table calculation
    #      detable = the energy step size used for array storage   
    #      emin = starting incident energy for table
    #      emax = mazximum incident energy for table calculation
    #	EP is ignored when iopt = 5 or 6
    # ************************************************************

    #      IN THE FOLLOWING THE LISTED VARIABLES ARE INDEXED
    #      THE INDICES I AND J STAND FOR THE FOLLOWING:
    #      - I -  SERIAL NUMBER OF ABSORBER LAYER (<20)
    #      - J -  SERIAL NUMBER OF ELEMENT WITHIN LAYER (<5)

    # *************************************************************


    #       ISG(I) = 0  - FOR SOLID ABSORBERS
    #          (I) = 1  - FOR GASEOUS ABSORBERS
    #       INN(I) = NUMBER OF ELEMENTS IN ONE LAYER
    #             (E.G. CH4 HAS TWO ELEMENTS # AND H)
    #       DEN(I) = DENSITY OF ABSORBE (FOR SOLIDS)
    #       THK(I) = THICKNESS OF ABSORBER IN MG/CMSQ (FOR SOLIDS)
    #       PRS(I) = PRESSURE (IN MM HG) FOR GAS ABSORBER
    #       XLN(I) = PHYSICAL LENGTH OF ABSORBER ( IN CM ) FOR GAS

    #  **** ABSORBER COMPOSITION ****

    #      ELNUM(I,J) = NUMBER OF ATOMS OF ELEMENT J IN LAYER I
    #      CONCN(I,J) = CONCENTRATION (MOLAR) OF ELEMENT J IN LAYER I
    #      ANUMB(I,J) = MASS   NUMBER   OF ELEMENT J IN LAYER I
    #      ZNUMB(I,J) = ATOMIC NUMBER   OF ELEMENT J IN LAYER I

    #      PRDEN(I,J) = PARTIAL DENSITY OF ELEMENT J IN LAYER I
    #      PRTHK(I,J) = PARTIAL THICKNESS OF ELEMENT J IN LAYER I (MG/CMSQ)

    #      MAXIMUM OF NINETEEN LAYERS SPECIFIED


    # *************************************************************
    # Loop to read in absorber parameters from a file...mostly
    for i = 1:ianz
        if ISG[i] != 1
            TOUT[i] = THK[i] / (DEN[i] * 1000)
            TOUTE[i] = TOUT[i] / 2.54
            break
        end
    end

    for i = 1:ianz
        INNW = INN[i]
        for j = 1:INNW
            ANUMBW[j] = ANUMB[i, j]
            ZNUMBW[j] = ZNUMB[i, j]
            ELNUMW[j] = ELNUM[i, j]
            CONCNW[j] = CONCN[j]
        end
        DENW = DEN[i]
        XLNW = XLN[i]
        PRSW = PRS[i]
        THKW = THK[i]
        if ISG[i] != 1
            #setabs(INNW, ANUMBW, ZNUMBW, ELNUMW, THKW, DENW)
            setlayer!(sandwich, ANUMBW, ZNUMBW, ELNUMW, THKW, DENW)
            for j = 1:INNW
                PRDEN[i, j] = PRDENW[j]
                PRTHK[i, j] = PRTHKW[j]
            end
        else
            setabg(INNW, ANUMBW, ZNUMBW, ELNUMW, CONCNW, THKW, DENW, PRSW, XLNW)
            DEN[i] = 0.0
            THK[i] = 0.0
            for j = 1:INNW
                PRDEN[i, j] = PRDENW[j]
                PRTHK[i, j] = PRTHKW[j]
                DEN[i] += PRDEN[i, j]
                THK[i] += PRTHK[i, j]
            end
        end
    end

    # *************************************************************

    #      START CALCULATION AND DETAILED PRINTOUT

    # *************************************************************
    #
    if iopt >= 5
        ep = emintab
        zp = zth + 1.0
        indexz = Int(zp - zth + 0.001)
    end
    #
    #299	continue   !  come here for new particle (zp change)

    izp = Int(zp + 0.001)
    #
    if iopt >= 5
        if izp > 70
            # write (6,*) 'no mass for Z = ',izp
            # stop
        else
            ap = amass(izp)

            # The trick!. To calculate energy losses for deuterons and tritons,
            # enter zth = 0.2 and 0.3 respectively.    E. Chavez jul/92

            if izp == 1
                iap = Int(zth * 10.0 + 0.1)
                ap = float(iap)
            end
        end
    end
    #
    if iopt == 6
        ideltalay = 8  # choose DE3 for eloss signal
        if izp == 1
            ideltalay = 16  # choose DEh for eloss signal
        end
    end
    #
    #300	CONTINUE    ! come here for new energy (ep)
    #
    EI = ep
    XUPDN = -1.0
    EPS = 0.0001
    I1STPASS = 1

    if iopt == 4
        #GO TO 600
    end

    ipunch = 2
    for i = 1:IANZ
        if iopt >= 5
            # GO TO 504
            # 504 skips the writing just below
        end

        # XNS - initial no. of intervals for integration of DE
        XNS = 2
        # EI = energy in
        ads(I, XUPDN, XNS, EPS, ap, zp, EI, DEI, ISTAT)
        # DEI = energy out - energy in (< 0.0 for energy loss)
        EIOLD = EI
        # E(I) = energy left after I'th element (EP-DE(1)-DE(2)+...)
        # if particle stopped in detector this is equal to energy lost
        # in remaining layers
        DE[i] = DEI
        E[i] = EI + DEI
        EI = E[i]
        INS = Int(XNS + 0.001)

        if iopt >= 5
            # GO TO 505
            # 505 skips loste and if statement below
        end
        # WRITE(IO2,401) INS,EIOLD,EI
        loste[i] = -1 * de[i]
        if EI < EPS || ISTAT == -1
            # WRITE(IO2,402) I
        end
        # 505
        loste[ianz+1] = e[ianz]
        # if particle stopped in layer beyond ianzi we must 
        # check if iopt=5 or6 and calculate the energy loss in the
        # front part (layers 1 thru ianzi).
        istore = i
        if EI < EPS
            ipunch = 1
        end
        # control loop exit
        # exit when particle runs out of energy in last layer
        if EI < EPS || ISTAT == -1
            # go to 701
        end
        # if interested in DE signal only (iopt=6) exit after 
        # layer for which DE is sought was traversed
        if iopt == 6 && I >= ianzide
            # go to 701
        end
    end

    #	this part for iopt=5,6 stores incident energy values

    # 701	continue
    if iopt != 5 && iopt != 6
        # go to 520
    end

    #  establish higher energy cutoof for next step with higher z
    #  should save time in calulating energies of particles stopped
    #  in the dead layer
    if I <= ianzi
        emintabz[izp] = ep
    end

    #  evalue = energy of particle when entering sensitive volume of det.
    #           when particle is stopped (ipunch=1) this is what is left
    #           once you take off the energy lost in the dead layer.
    evalue = E[ianzi]
    #         = when particle punches thouough detector its energy at the
    #           end is EI - subtract this from what it entered with (after
    #           dead layers) and again you got the energy deposited.
    if ipunch == 2
        evalue = E[ianzi] - EI >= 0.0 ? E[ianzi] - EI : 0.0
        #      roundoff errors could resutl in negative energies !!
    end

    #  EI is current energy after last layer - is nonzero when particle
    #  punched through and to get energy deposited must subtract this 
    #  "left over" energy from the energy of particle had when it entered
    #  the detector's sensitive volume
    indexe = Int((ep + 0.001) / detable) + 1
    indexz = Int(zp - zth + 0.001)
    if iopt == 5
        eptable[indexz, indexe, ipunch] = evalue
    end
    if iopt == 6
        ipunch = 1
        eptable[indexz, indexe, ipunch] = -DE[ianzide]
    end

    #	now repeat calculation for same Z but new energy
    ep = ep + detable
    if ep > emaxtab
        # go to 709
    end
    # go to 300

    # 709	continue
    #	reset energy to emintab and up zp by one until we top
    #	zmax - this portion controls looping over Z!
    for ieps = 1:ianz
        E[ieps] = 0.0
        DE[ieps] = 0.0
    end

    ep = emintabz[izp]
    zp = zp + 1.0
    izp = Int(zp + 0.001)
    emintabz[izp] = emintabz[izp-1]

    eltime = eltimeo #secnds (eltimeo)
    itminutes = Int(eltime / 60.0)
    tminutes = float(itminutes)
    tseconds = eltime - tminutes * 60.0
    eltimeo = eltimeo + eltime
    zpm1 = zp - 1.0

    if izp >= izmax
        # go to 711
    end
    # go to 299

    # 711	continue

    #	get here when iopt = 5 or 6 calculation is done
    #	now ready to store array on disk
    # write(io0,712) zth,zmax,emintab,emaxtab,detable

    iemintab = Int((emintab + 0.001) / detable) + 1
    iemaxtab = Int((emaxtab + 0.001) / detable)
    iztop = Int(zp - 1.0 - zth + 0.001)

    for ipp = 1:2
        for indexz = 1:iztop
            iztab = indexz + Int(zth + 0.01)
            imasstab = Int(amass[iztab] + 0.1)
            if iztop == 1
                imasstab = Int(ap + 0.01)
            end
            # write(io0,7712) indexz,iztab,imasstab,iemintab,iemaxtab
            # 7712	format(5i20)
            itblow = iemintab
            for indexe = iemintab:iemaxtab:10
                itbup = itblow + 9
                # write(io0,713) (eptable(indexz,itbe,ipp),itbe=itblow,itbup)
                itblow = itbup + 1
            end
        end
    end

    #      MORE INPUT FOR NEW CALCULATION WITH SAME ABSORBERS

    # 520	CONTINUE


    zp = -1.0
    if zp <= 0.0
        # GO TO 2000
    end
    izp = Int(zp + 0.001)
    if ap <= 0.0
        AP = amass[izp]
    end
    # GO TO 299

    # 600
    for i = 1:ianz
        ILAY = i
        XNS = 2.0
        ads(i, XUPDN, XNS, EPS, ap, zp, EI, DEI, ISTAT)
        EIOLD = EI
        DE[i] = DEI
        E[i] = EI + DEI
        # E(I) = energy left after I'th element (EP+DE(1)+DE(2)+...)
        # if particle stopped in detector this is equal to energy lost
        # in remaining layers
        xmem[i] = E[i]
        EI = E[i]
        INS = Int(XNS + 0.001)
        if EI <= 0.0 || ISTAT == -1
            # GO TO 601
        end
    end

    # 601
    if ISTAT == 0
        if EI < 0.003 && EI >= 0.0
            izp = Int(zp + 0.001)
            if ap <= 0.0
                ap = amass(izp)
            end
            if zp >= 0.0 # else go to 2000
                IPASS = 0
                I1STPASS = 1
                EI = ep
                # GO TO 600
            end
        end
    else
        isign = -isold
        isold = -1
    end

    if I1STPASS > 0
        isign = 1
        I1STPASS = 0
        if ISTAT == 0
            DSEP = -DSEP0
        else
            DSEP = DSEP0
        end
    end

    #	IF THE INITIAL ENERGY WAS TOO LARGE, THEN THE ION WILL PUNCH THROUGH
    #	THE DETECTOR A NUMBER OF TIMES UNTIL THE ENERGY IS REDUCED BELOW 
    #	THE PUNCH-THROUGH ENERGY (PTE). IF THE INITIAL ENERGY WAS TOO SHORT
    #	THEN IT WON'T UNTIL PTE IS REACHED. IN THIS MOMENT, IPASS IS SET TO
    #	ONE, AND FROM THIS POINT EVERY FURTHER CALCULATION WILL IMPLY A
    #	REDUCTION BY HALF OF THE SIZE OF "DSEP", AND A CHANGE OF SIGN ONLY
    #	IF APROPRIATE.

    if isign < 0
        IPASS = 1#IPASS=1 UNTIL "PTE" IS FOUND.
    end
    if IPASS == 1
        DSEP = isign * DSEP * 0.5
    end

    if abs(DSEP) < 0.05
        EI = 0.00001
        ISTAT = 0
        # GO TO 601
    end

    EP = EP + DSEP
    ei = ep

    # GO TO 600

    # 1000	CONTINUE
    # GO TO 10
    # 2000	CONTINUE
    return
end
=#

function stoppingenergy(A::Integer, Z::Integer, sandwich::Sandwich;
    initialenergy::Unitful.Energy=100.0u"MeV", integrationsteps::Integer=1000,
    precision::Float64=1e-4, atol::Real=1e-3)
    energy = uconvert(1.0u"MeV", initialenergy)
    atol *= 1.0u"MeV"
    stopped = false
    Eloss = 0.0u"MeV"

    while !stopped
        stopee = Particle(A, Z, energy)
        try
            @inbounds for layer in sandwich.layers
                Eloss += ads(stopee, layer, integrationsteps, precision)
            end
            ΔE = energy - Eloss
            if isapprox(ΔE, 0.0u"MeV"; atol=atol)
                return energy
            end
            # Reset variables and lower energy
            energy -= ΔE
            Eloss = 0.0u"MeV"
        catch err
            if isa(err, ParticleStoppedException)
                # If the particle stopped need to increase the energy and try again
                energy *= 2
                Eloss = 0.0u"MeV"
            else
                rethrow()
            end
        end
    end
end

function stoppingdepth(stopee::Particle, layer::T; integrationsteps::Integer=1000, precision::Float64=1e-4) where {T<:AbstractAbsorber}

end

function dedx(energy::Unitful.Energy, index::Integer, layer::AbstractAbsorber, part::Particle)
    δEδx = 0.0u"MeV*cm^2/mg"

    vel = _velocity(part)
    XI = vel^2 / layer.Z[index]
    AZ1 = log(1.035 - 0.4 * exp(-0.16 * part.Z))
    VZ1 = (-116.79 - 3350.4 * vel * part.Z^(-0.509)) * vel * part.Z^(-0.509)
    DEDXHIfactor = part.Z
    if VZ1 > -85.2
        DEDXHIfactor = part.Z * (1 - exp(VZ1))
    end
    if part.Z > 2
        FV = vel <= 0.62 ? 1 - exp(-137.0 * vel) : 1.0
        DEDXHIfactor = part.Z * (1 - exp(FV * AZ1 - 0.879 * 137.0 * vel / (part.Z^(0.65))))
    end

    # ENERGY LOSS OF HEAVY IONS
    # EFFECTIVE CHARGE
    # ELECTRONIC ENERGY LOSS DEDXHI
    DEDXHI = DEDXHIfactor^2 * layer.Z[index] * _Y(XI, index, layer) /
             (1.0u"u^-1" * layer.mass[index] * vel^2)

    # NUCLEAR ENERGY LOSS DEDXNU
    ZA = sqrt(part.Z^(2 / 3) + layer.Z[index]^(2 / 3))

    ϵ = 3.25e4u"MeV^-1" * layer.mass[index] * energy /
        (part.Z * layer.Z[index] * (part.mass + layer.mass[index]) * ZA)

    Σ = 1.7 * sqrt(ϵ) * log(ϵ + 2.1718282) / (1 + 6.8 * ϵ + 3.4 * ϵ^1.5)

    DEDXNU = Σ * 5.105u"u" * part.Z * layer.Z[index] * part.mass /
             (ZA * layer.mass[index] * (part.mass + layer.mass[index]))

    # TOTAL ENERGY LOSS
    δEδx += (DEDXHI + DEDXNU) * 1.0u"MeV*cm^2/mg"
    return δEδx
end

function _Y(XI::Float64, index::Integer, layer::SolidAbsorber)
    pdens = layer.partialdensity[index] * 1.0u"cm^3/g"
    FY = 54721.0 * (1 + 0.0515u"u^(-1/2)" * sqrt(layer.mass[index] / pdens) - exp(-0.23 * layer.Z[index]))

    # CALCULATION OF Y(XI)
    Y = 3.3e-4 * log(1 + XI * FY)

    if 1e-9 <= XI <= 5e-4
        FG = 1.2E-4 * layer.Z[index]^2 + 2.49E-2u"u^-1" * layer.mass[index] / pdens
        HZ2 = 1.32e-5 * (9.0 - (_G1(layer.Z[index]) + _G2(layer.Z[index]) + _G3(layer.Z[index]) + _G4(layer.Z[index]) + _G5(layer.Z[index])))
        C = 2 * sqrt(XI) / (layer.Z[index] * (1 + 1.E4 * sqrt(XI)))
        AL = log(XI * FG / 2.7E-5)
        Y += (C - HZ2 * AL * exp(-0.32 * AL^2)) / (1 + (XI * 1E4)^3)
    end
    return Y
end

function _Y(XI::Float64, index::Integer, layer::GasAbsorber)
    FY = 54721.0 * (1.35 - exp(-0.13 + 0.0014 * layer.Z[index]))

    # CALCULATION OF Y(XI)
    Y = 3.3E-4 * log(1 + XI * FY)

    if 1e-9 <= XI <= 5e-4
        FG = 1.3 / (1 + exp(3.0 - layer.Z[index] / 5.0))
        HZ2 = 1.32e-5 * (9.0 - (_G1(layer.Z[index]) + _G2(layer.Z[index]) + _G3(layer.Z[index]) + _G4(layer.Z[index]) + _G5(layer.Z[index])))
        C = sqrt(XI) / (layer.Z[index] * (1 + 1.E4 * sqrt(XI)))
        AL = log(XI * FG / 2.7E-5)
        Y += (C - HZ2 * AL * exp(-0.32 * AL^2)) / (1 + (XI * 1E4)^3)
    end
    return Y
end

@inline function _G1(Z2::Integer)
    if Z2 <= 26
        return 19.84 * exp(-0.17 * (Z2 - 4.25)^2)
    end
    return 0.000001
end

@inline function _G2(Z2::Integer)
    if Z2 <= 38
        return 17.12 * exp(-0.12 * (Z2 - 11.63)^2)
    end
    return 0.000001
end

@inline function _G3(Z2::Integer)
    return 7.95 * exp(-0.015 * (Z2 - 30.2)^2)
end

@inline function _G4(Z2::Integer)
    5.84 * exp(-0.022 * (Z2 - 48.63)^2)
end

@inline function _G5(Z2::Integer)
    7.27 * exp(-0.005 * (Z2 - 73.06)^2)
end

function ads(part::Particle, layer::T, steps::Integer, precision::Float64) where {T<:AbstractAbsorber}
    Eᵢ = part.energy  # Energy at step i
    δEₜ = 0.0u"MeV"  # Total energy loss
    δEᵢ = 0.0u"MeV"  # Energy loss from first step
    δEₙ = 0.0u"MeV"

    i = 1
    while i < steps
        try
            δEₙ = _dedxloop!(Eᵢ, i, steps, layer, part)
        catch err
            if isa(err, ZeroEnergyException)
                i = 1
                steps *= 2
                Eᵢ = part.energy
                δEₜ = 0.0u"MeV"
            else
                rethrow()
            end
        else
            if i == 1
                δEᵢ = δEₙ
            elseif i == 2
                calcprecision = abs(δEᵢ - δEₙ) / (δEᵢ + δEₙ)
                if calcprecision > precision
                    i = 1
                    steps *= ceil(Int, calcprecision / precision)
                    Eᵢ = part.energy
                    δEₜ = 0.0u"MeV"
                    continue
                end
            end
            δEₜ += δEₙ
            Eᵢ -= δEₙ
            i += 1
        end
    end

    return δEₜ
end

function _dedxloop!(Eᵢ::typeof(1.0u"MeV"), index::Integer, steps::Integer, layer::T, part::Particle) where {T<:AbstractAbsorber}
    δE = 0.0u"MeV"
    @inbounds for (j, pthick) in enumerate(layer.partialthickness)
        thickstep = pthick / steps
        δEδx = dedx(Eᵢ, j, layer, part)
        if Eᵢ - δEδx * thickstep < 0.0u"MeV"
            if index <= 2
                steps *= 2
                index = 1
                Eᵢ = part.energy
                throw(ZeroEnergyException())
            else
                throw(ParticleStoppedException())
            end
        end
        δE += δEδx * thickstep
    end
    return δE
end

@inline function _velocity(mass::typeof(1.0u"u"), energy::typeof(1.0u"MeV"))
    return sqrt(2.13E-3u"u/MeV" * energy / mass)
end

@inline function _velocity(part::Particle)
    return _velocity(part.mass, part.energy)
end

function fkinem(EP, AP, AT, TH)
    if AP <= AT
        E = EP * AP^2 / (AP + AT)^2
        E = E * (cos(TH) + sqrt((AT / AP)^2 - sin(TH)^2))^2
        return E
    end
    return 0.0
end

function ffkin(EP, AP, AT, TH, Q)
    # 	INELASTIC SCATTERING
    B = AP^2 * EP / (AP + AT)^2 / (EP + Q)
    D = AT^2 / (AP + AT)^2 * (1 + AP * Q / AT / (EP + Q))
    E = (EP + Q) * B * (COS(TH) + SQRT(D / B - SIN(TH)^2))^2
    return E
end

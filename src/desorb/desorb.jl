export desorb!

include("absorbers.jl")
include("standardmedia.jl")

function desorb!(ianz, zp, ap, ep, loste)

    # **  CALCULATES ENERGY LOSS IN AN ABSORBER SANDWICH *******
    # ****  ENERGY DEPOSIT IN SECTIONS OF IONIZATION  **********
    # ***************  CHAMBER (DE1, DE2, AND DE3)  ************


    #	PARAMETER (DSEP0=4.0,isold0=1)

    #       COMMON ISG(19),INN(19),DEN(19),THK(19),PRS(19),XLN(19),ARDEN(19)
    #	1      ,ZNUMB(19,4),ANUMB(19,4),ELNUM(19,4),CONCN(19,4)
    #	1         ,PRDEN(19,4),PRTHK(19,4)
    #	DIMENSION ZNUMBW(4),ANUMBW(4),ELNUMW(4),CONCNW(4)
    #	1         ,PRDENW(4),PRTHKW(4)
    #	DIMENSION E(19),DE(19),xmem(19)
    #	DIMENSION TOUT(19),TOUTE(19)
    #	dimension eptable(50,500,2),emintabz(50),ianzv(5)
    #	real loste(19)
    #	DATA io1,IO2,IO3,io0/9,11,12,1/
    #	data iopt,ianzi,ianzide/1,2,1/
    #
    #	save eltimo,izmax,izth

    #       the mass table is to be used only for iopt = 5,6
    #       use atomic masses to average for isotipic composition.
    #       taken from Formulas Facts and Constants, H. J. Fischbeck and
    #       K. H. Fischbeck. Springer - Verlag 1987 2nd ed, pages 164-183.

    #        dimension amass(70)
    #        data amass/1.01,4.00,6.94,9.01,10.81,12.01,14.01,16.00,19.00
    #     1,20.18,22.99,24.31,26.98,28.09,30.97,32.07,35.45,39.95
    #     2,39.10,40.08,44.96,47.88,50.94,52.00,54.94,55.85,58.93
    #     3,58.69,63.55,65.39,69.72,72.59,74.92,78.96,79.90,83.80
    #     4,85.47,87.62,88.91,91.22,92.91,95.94,98.,101.07,102.91
    #     5,106.42,107.87,112.41,114.82,118.71,121.75,127.60,126.90
    #     6,131.29,132.91,137.33,138.91,140.12,140.91,144.24,147.
    #     7,150.36,151.96,157.25,158.93,162.5,164.93,167.26,168.93
    #     8,173.04/

    #       for Z > 70, you are in trouble !!
    #
    #      open(unit=IO1, status='OLD', file='absorb.inp')
    #      open(unit=IO2, status='OLD', file='desorb.out')
    #      rewind IO1
    #      rewind IO2
    #
    #10	CONTINUE


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

    #	read(io1,*) iopt
    #c1	FORMAT(10I)
    #	DSEP = DSEP0
    #	isold = isold0
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
    #	read(IO1,*) ianz,ianzi,ianzide
    #
    #	ianzv(5)
    #	     = Index of layer where exiting particle velocit is 
    #	       calculated (only for option 3) max = 5
    #
    #	if(iopt.eq.3) read (IO1,*) knz
    #	if(iopt.eq.3) read (IO1,*) (ianzv(k),k=1,knz)
    # ***********************************************************


    #      AP = PROJECTILE MASS
    #      ZP = PROJECTILE CHARGE
    #      EP = PROJECTILE ENERGY
    #
    #	read(io1,*) zp,ap,ep
    #2	FORMAT(f10.3)

    # *************************************************************
    #
    #      zth= threshold Z for incident energy tble calc. (iopt=5,6)
    #      zmax = maximum z for table calculation
    #      detable = the energy step size used for array storage   
    #      emin = starting incident energy for table
    #      emax = mazximum incident energy for table calculation
    #	EP is ignored when iopt = 5 or 6
    #
    #	if(iopt.ge.5) then
    #		read(IO1,*) zth,zmax,emintab,emaxtab,detable
    #		eltimeo = secnds(0.0)  !  start timing
    #		izth = ifix (zth + 0.01)
    #		izmax = ifix (zmax + 0.01)
    #
    #	if detable > emintab change emintab to 1/2 * detable
    #
    #		if(emintab.lt.detable) emintab = 0.5*detable
    #		write(io2,291) zth+1,zmax,emintab,emaxtab,detable
    #	else
    #	endif
    #
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
        for j = 1:INN[i]
            # read(io1,*) ANUMB(ISN,IMN),ZNUMB(ISN,IMN),ELNUM(ISN,IMN),CONCN(ISN,IMN)
        end
    end
    #	DO 100 ISN=1,IANZ
    #		read(io1,*) ISG(ISN),INN(ISN)
    #		IF(ISG(ISN).EQ.1) GO TO 50
    #		read(io1,*) DEN(ISN),THK(ISN)
    #		TOUT(ISN)=THK(ISN)/(DEN(ISN)*1000.)
    #		TOUTE(ISN)=TOUT(ISN)/2.54
    #		GO TO 100
    #50              continue
    #		read(io1,*) PRS(ISN),XLN(ISN)
    #		DO 60 IMN=1,INN(ISN)
    #			read(io1,2) ANUMB(ISN,IMN),ZNUMB(ISN,IMN),
    #	1		ELNUM(ISN,IMN),CONCN(ISN,IMN)
    #60		CONTINUE
    #100	CONTINUE

    # ****************************************************************

    #	if(iopt.ne.5.and.iopt.ne.6) WRITE(IO2,101) ap,zp,ep
    #	if(iopt.ne.5.and.iopt.ne.6) WRITE(IO2,102) ianz

    # *****************************************************************
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

            #       The trick!. To calculate energy losses for deuterons and tritons,
            #       enter zth = 0.2 and 0.3 respectively.    E. Chavez jul/92

            if izp == 1
                iap = Int(zth * 10.0 + 0.1)
                ap = float(iap)
            end
        end
    end
    #
    if iopt == 6
        ideltalay = 8                # choose DE3 for eloss signal
        if izp == 1
            ideltalay = 16   # choose DEh for eloss signal
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
        if ISG[i] == 0
            # WRITE(IO2,311) I,THK(I),TOUT(I),TOUTE(I),DEN(I)
        elseif ISG[i] == 1
            # WRITE(IO2,312) I,THK(I),PRS(I),XLN(I),DEN(I)
        end
        for j = 1:INN[i]
            # WRITE(IO2,321) ANUMB(I,J),ZNUMB(I,J),PRTHK(I,J)
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
    if (izp - 1) <= izmax
        # write(io2,714) zpm1,ap,tminutes,tseconds
        # 714 format(' finished Z=',f3.0,'   A=',f4.0,' -  ',f5.0,
        # +' minutes and ',f3.0,' seconds elapsed')
    end

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
    # 712	format(5e16.8)
    # 713	format(10e16.8)
    # close (unit = io0)

    #      MORE INPUT FOR NEW CALCULATION WITH SAME ABSORBERS

    # 520	CONTINUE


    #	read(io1,*) zp,ap,ep
    zp = -1.0
    if zp <= 0.0
        # GO TO 2000
    end
    izp = Int(zp + 0.001)
    if ap <= 0.0
        AP = amass[izp]
    end
    if zp > 0.0
        #WRITE(IO2,101) ap,zp,ep
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
        # WRITE (IO2,613) ILAY,INS,EI,ISTAT
        if EI <= 0.0 || ISTAT == -1
            # GO TO 601
        end
    end
    #600	DO I = 1, ianz
    #		ILAY = I
    #		XNS = 2.0
    #		CALL ADS(I,XUPDN,XNS,EPS,ap,zp,EI,DEI,ISTAT)
    #		EIOLD = EI
    #		DE(I) = DEI
    #		E(I) = EI + DEI
    #       E(I) = energy left after I'th element (EP+DE(1)+DE(2)+...)
    #       if particle stopped in detector this is equal to energy lost
    #       in remaining layers
    #		xmem(i) = E(I)
    #		EI = E(I)
    #		INS = IFIX(XNS + 0.001)
    #		WRITE (IO2,613) ILAY,INS,EI,ISTAT
    #		IF (EI.LE.0.0.OR.ISTAT.EQ.-1) GO TO 601
    #	END DO

    # 601
    if ISTAT == 0
        if EI < 0.003 && EI >= 0.0
            # WRITE(IO2,611) ap,zp,ep,ILAY,xmem(5),xmem(6),xmem(7)
            # read(io1,*) zp,ap
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

    #601	IF (ISTAT.EQ.0)THEN
    #		IF (EI.LT.0.003.AND.EI.GE.0.0) THEN
    #	WRITE(IO2,611) ap,zp,ep,ILAY,xmem(5),xmem(6),xmem(7)
    #		read(io1,*) zp,ap
    #		izp = ifix (zp + 0.001)
    #		if(ap.le.0.) ap = amass(izp)
    #			IF (zp.LT.0.0) GO TO 2000
    #			IPASS = 0
    #			I1STPASS = 1
    #			EI = ep
    #			GO TO 600
    #		END IF
    #		isign = isold
    #		isold = 1
    #	ELSE
    #		isign = - isold
    #		isold = -1
    #	END IF

    #	IF (I1STPASS.gt.0) THEN
    #		isign = 1
    #		I1STPASS = 0
    #		IF (ISTAT.EQ.0) THEN
    #			DSEP =  - DSEP0
    #		ELSE
    #			DSEP = DSEP0
    #		END IF
    #	END IF

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

    #101	FORMAT(//////' PASSAGE OF CHARGED PARTICLE THROUGH ABSORBER',
    #	1       ' SANDWICH '////'   AP = ',F6.0,'   ZP = ',F5.0,
    #	1       '    INITIAL ENERGY = ',F12.5//)
    #102	FORMAT('   ABSORBER SANDWICH CONTAINS - ',I2,' LAYERS'//)
    #291	format(' start absorber calculations for   z =',f3.0
    #     +,'  to z =',f3.0,/' and energies from emin =',f6.2
    #     +,'  to emax =',f7.2,'  in ',f6.2,'MeV steps')
    #311	FORMAT(//'  LAYER # ',I2,'    - SOLID ABSORBER -   ',
    #	1' AREAL DENSITY = ',E10.4/' THICKNESS = ',E10.4,
    #	1' CM   OR ',E10.4,' INCH      DENSITY =',E10.3,' G/CM3')
    #312	FORMAT(//'   LAYER # ',I2,'    -  GAS  ABSORBER -   ',
    #	1'  AREAL DENSITY = ',E10.4/'  PRESSURE = ',E10.4,
    #	1' TORR    LENGTH ',E9.4,'CM    DENSITY =',E9.3,'MG/CM3')
    #321	FORMAT(7X,' A =',F6.0,'  Z =',F5.0,'   AREAL DENSITY'
    #	1,' (PARTIAL) = ',E12.5,' MG/CMSQ')
    #401	FORMAT(' CALC IN-',I4,' STEPS'
    #	1,'   ENERGY IN = ',F8.3,'    ENERGY OUT = ',F8.3
    #	2,'(MEV)')
    #402	FORMAT(' CHARGED PARTICLE STOPPED IN LAYER # ',I2)
    #613	FORMAT (2X,'LAYER= ',I2,': ',I6,' ITERATIONS'
    #	1,   ', E final= ',F10.4,' STATUS= ',I2)
    #611	FORMAT (2X,'Ion  (A , Z): (',F4.0,' , ',F3.0
    #	1,'), E(MeV)= ',F7.2,'  STOPPED IN LAYER ',I2/ 
    #	1'   Esum = ',f7.2,'    Esum-E1 = ',f7.2,
    #	1'   Esum-E1-E2 = ',f7.2)
    #703	format('  eptable (',i2,', ',i3,', ',i1,' ) =',f8.2)'
    return
end

function dedx(layer::T, stopee::Particle) where {T<:AbstractAbsorber}
    DEDXTO = 0.0
    vel = velocity(stopee)
    for i in 1:length(layer.A)
        XI = vel^2 / 2

        # ENERGY LOSS OF HEAVY IONS
        # EFFECTIVE CHARGE
        # ELECTRONIC ENERGY LOSS DEDXHI
        DEDXHI = layer.Z[i] * _Y(XI, i, layer) / (A2 * vel^2)

        AZ1 = log(1.035 - 0.4 * exp(-0.16 * stopee.Z))

        if stopee.Z > 2
            FV = vel < 0.62 ? 1 - exp(-137.0 * vel) : 1.0
            DEDXHI *= (stopee.Z * (1 - exp(FV * AZ1 - 0.879 * VV0 / (Z1^(0.65)))))^2
        elseif VZ1 > -85.2
            DEDXHI *= (stopee.Z * (1 - exp(VZ1)))^2
        else
            DEDXHI *= stopee.Z^2
        end

        # NUCLEAR ENERGY LOSS DEDXNU
        ZA = sqrt(stopee.Z^(2 / 3) + layer.Z[i]^(2 / 3))
        ϵ = 3.25E4 * layer.A[i] * stopee.energy / (stopee.Z * layer.Z[i] * (stopee.A + layer.A[i]) * ZA)
        Σ = 1.7 * sqrt(ϵ) * log(ϵ + 2.1718282) / (1 + 6.8 * ϵ + 3.4 * ϵ^1.5)
        DEDXNU = Σ * 5.105 * stopee.Z * layer.Z[i] * stopee.A / (ZA * layer.A[i] * (stopee.A + layer.A[i]))

        # TOTAL ENERGY LOSS DEDXTO

        DEDXTO += DEDXHI + DEDXNU
    end
    return DEDXTO
end

function _Y(XI::Float64, index::Integer, layer::SolidAbsorber)
    FY = 54721.0 * (1.35 - exp(-0.13 + 0.0014 * layer.Z[index]))
    FG = 1.3 / (1 + exp(3.0 - layer.Z[index] / 5.0))

    HZ2 = 1.32e-5 * (9.0 - (_G1(layer.Z[index]) + _G2(layer.Z[index]) + _G3(layer.Z[index]) + _G4(layer.Z[index]) + _G5(layer.Z[index])))

    # CALCULATION OF Y(XI)
    Y = 3.3E-4 * log(1 + XI * FY)

    if 1.E-9 <= XI <= 5.E-4
        C = 2 * sqrt(XI) / (layer.Z[i] * (1 + 1.E4 * sqrt(XI)))
        AL = log(XI * FG / 2.7E-5)
        Y += (C - HZ2 * AL * exp(-0.32 * AL^2)) / (1 + (XI * 1E4)^3)
    end
end

function _Y(XI::Float64, index::Integer, layer::GasAbsorber)
    FY = 54721.0 * (1 + 0.0515 * sqrt(layer.A[index] / layer.partialdensity[index]) - exp(-0.23 * Z2))
    FG = 1.2E-4 * layer.Z[index]^2 + 2.49E-2 * layer.A[index] / layer.partialdensity[index]

    HZ2 = 1.32e-5 * (9.0 - (_G1(layer.Z[index]) + _G2(layer.Z[index]) + _G3(layer.Z[index]) + _G4(layer.Z[index]) + _G5(layer.Z[index])))

    # CALCULATION OF Y(XI)
    Y = 3.3E-4 * log(1 + XI * FY)

    if 1.E-9 <= XI <= 5.E-4
        C = sqrt(XI) / (layer.Z[i] * (1 + 1.E4 * sqrt(XI)))
        AL = log(XI * FG / 2.7E-5)
        Y += (C - HZ2 * AL * exp(-0.32 * AL^2)) / (1 + (XI * 1E4)^3)
    end
end

function _G1(Z2::Integer)
    if Z2 <= 26
        return 19.84 * exp(-0.17 * (Z2 - 4.25)^2)
    end
    return 0.000001
end

function _G2(Z2::Integer)
    if Z2 <= 38
        return 17.12 * exp(-0.12 * (Z2 - 11.63)^2)
    end
    return 0.000001
end

function _G3(Z2::Integer)
    return 7.95 * exp(-0.015 * (Z2 - 30.2)^2)
end

function _G4(Z2::Integer)
    5.84 * exp(-0.022 * (Z2 - 48.63)^2)
end

function _G5(Z2::Integer)
    7.27 * exp(-0.005 * (Z2 - 73.06)^2)
end

function ads(I1, SIGN, XN1, EPS, A, Z, E, ISTAT)
    #	SUBROUTINE FOR ENERGY LOSS CALCULATIONS
    #	CALL DEDX FRO STOPPING POWER CALCULATIONS

    #	COMMON ISG(19),INN(19)
    #	1      ,DEN(19),THK(19),PRS(19),XLN(19),ARDEN(19)
    #	1      ,ZNUMB(19,4),ANUMB(19,4),ELNUM(19,4),CONCN(19,4)
    #	1         ,PRDEN(19,4),PRTHK(19,4)

    #	N1= NUMBER OF SUBDIVISIONS FOR INTEGRATION
    #		OF ENERGY LOSS

    #1000	CONTINUE
    EH = E
    N1 = Int(XN1 + 0.001)
    DEDNEXT = 0.0
    for i = 1:N1
        J1 = INN[I1]
        ISGW = ISG[I1]
        for j = 1:J1
            AX = ANUMB[I1, j]
            ZX = ZNUMB[I1, j]
            FX = PRTHK[I1, j]
            DENST = PRDEN[I1, j]
            VH = vel(EH, A)
            dedx(Z, A, ZX, AX, DENST, EH, VH, ISGW)
            EH += DE * SIGN * FX
            if EH <= 0.0
                if i <= 2
                    N1 *= 2
                    XN1 = float(N1)
                    # GO TO 1000
                else
                    ISTAT = -1
                    return
                end
            end
            if i <= 2
                DEDNEXT += DE * FX
            end
        end
        if i == 1
            DED1ST = DEDNEXT
            DEDNEXT = 0.0
        elseif i == 2
            DDD = DED1ST - DEDNEXT
            if DDD < 0
                DDD *= -1
            end
            DDS = DED1ST + DEDNEXT
            DDR = DDD / DDS
            if DDR > EPS
                N1 = N1 * 2
                XN1 = float(N1)
                # GO TO 1000
            end
        end
    end
    ISTAT = 0
    DEE = EH - E
    return DEE
end

function velocity(part::Particle)
    return sqrt(2.13E-3 * part.energy / part.A)
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

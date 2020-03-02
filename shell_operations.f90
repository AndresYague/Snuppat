MODULE shell_operations
    USE readvars_mod
    USE structures_mod
    USE math_routines
    IMPLICIT NONE
    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine inputs the chemical abundances of the first array into   !!!
!!! the second, interpolating where necessary.                               !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -chemArr, array with all the abundances on the star.                     !!!
!!! -model, array with all the physics in which to store the chemistry.      !!!
!!! -chemShells, size of chemArr.                                            !!!
!!! -siz, number of chemical species.                                        !!!
!!!                                                                          !!!
!!! On output, model has the interpolated abundance values.                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE chemInPhys(chemArr, model, chemShells, siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (SHELL), POINTER::chemArr(:), model(:)
    INTEGER::chemShells, siz
    
    ! Local
    DOUBLE PRECISION, ALLOCATABLE::mass(:), valMass(:), valDens(:, :), dens(:)
    INTEGER::ii, jj, modShells
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Prepare a mesh with the intermediate masses for both models and values for
    ! intermediate points. Then interpolate the new values and, finally, put
    ! them into the new model.
    modShells = SIZE(model)
    
    ! Allocate
    ALLOCATE(mass(modShells - 1), valMass(chemShells - 1))
    ALLOCATE(valDens(chemShells - 1, siz), dens(modShells - 1))
    
    ! Fill arrays
    DO ii = 1, modShells - 1
        mass(ii) = (model(ii)%mass + model(ii + 1)%mass)*0.5
    END DO
    DO ii = 1, chemShells - 1
        valMass(ii) = (chemArr(ii)%mass + chemArr(ii + 1)%mass)*0.5
        valDens(ii, :) = (chemArr(ii)%dens + chemArr(ii + 1)%dens)*0.5
    END DO
    
    ! Now perform interpolation for each species
    DO ii = 1, siz
        ! Interpolate
        ! The absolute value is used because it is possible to get a negative
        ! abundance when interpolating. It would be a spurious one, 20 orders
        ! of magnitude below those surrounding it.
        dens = ABS(interpolate(mass, valMass, valDens(:, ii), modShells - 1))
        
        ! Introduce values
        DO jj = 1, modShells
            IF (jj.EQ.1) THEN
                ! First case, we can put an arbitrary first value as all of the
                ! others will be based on this one and the real values will be
                ! the averages. We will choose the average value as its value.
                model(jj)%dens(ii) = dens(jj)
            ELSE
                ! In any other case we want the average equal to
                ! the interpolated value
                model(jj)%dens(ii) = 2*dens(jj - 1) - model(jj - 1)%dens(ii)
            END IF
        END DO
    END DO
    
    ! Deallocate
    DEALLOCATE(mass, valMass, valDens, dens)
    
END SUBROUTINE chemInPhys

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates and inputs the intershell values we are       !!!
!!! interested in, in intShell, returning also nIntShell, which is used.     !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -model, array with all the physics and chemistry.                        !!!
!!! -intShell, inter-shell array.                                            !!!
!!! -totShell, size of model.                                                !!!
!!! -siz, number of chemical species.                                        !!!
!!!                                                                          !!!
!!! On output, intShell has the intershell physics and chemistry.            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE storeShells(model, intShell, totShell, siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (SHELL)::model(:)
    TYPE (INTERSHELL), POINTER::intShell(:)
    INTEGER::totShell, siz
    
    ! Local
    DOUBLE PRECISION::avRad, dm
    INTEGER::ii
    LOGICAL::currentConv, nextConv
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Store all intershells
    ii = 1
    DO WHILE (ii.LT.totShell)
        ! Nullify convShell
        NULLIFY(intShell(ii)%convShell)
        
        ! Calculate average convective criterion
        avRad = model(ii + 1)%radiat + model(ii)%radiat
        
        ! Calculate mass extension
        dm = model(ii + 1)%mass - model(ii)%mass
        
        ! See if shell is convective (only if the mass is not zero)
        IF ((avRad.GT.0.D0).AND.(dm.GT.1.D-50)) THEN
            intShell(ii)%isConvective = .TRUE.
        END IF
        
        ! Calculate and store values
        intShell(ii)%temp = (model(ii + 1)%temp + model(ii)%temp)*0.5D0
        intShell(ii)%rho = (model(ii + 1)%rho + model(ii)%rho)*0.5D0
        
        intShell(ii)%rho0 = model(ii)%rho
        intShell(ii)%rho1 = model(ii + 1)%rho
        
        intShell(ii)%mass0 = model(ii)%mass
        intShell(ii)%mass1 = model(ii + 1)%mass
        
        intShell(ii)%rad0 = model(ii)%radius
        intShell(ii)%rad1 = model(ii + 1)%radius
        
        intShell(ii)%hp0 = model(ii)%hp
        intShell(ii)%hp1 = model(ii + 1)%hp
        
        intShell(ii)%pres0 = model(ii)%pressure
        intShell(ii)%pres1 = model(ii + 1)%pressure
        
        intShell(ii)%vel0 = model(ii)%velocity
        intShell(ii)%vel1 = model(ii + 1)%velocity
        
        intShell(ii)%radiat0 = model(ii)%radiat
        intShell(ii)%radiat1 = model(ii + 1)%radiat
        
        ! Store abundances
        ALLOCATE(intShell(ii)%dens(siz))
        intShell(ii)%dens = (model(ii + 1)%dens + model(ii)%dens)*0.5D0
        
        ii = ii + 1
    END DO
    
    ! Now check if there's a lone radiative or convective shell and jump over it
    DO ii = 1, totShell - 2
        currentConv = intShell(ii)%isConvective
        nextConv = intShell(ii + 1)%isConvective
        
        ! Check if first shell is a lone convective shell
        IF (ii.EQ.1) THEN
            IF (currentConv.AND.(.NOT.nextConv)) THEN
                IF (.NOT.intShell(ii + 2)%isConvective) THEN
                    intShell(ii)%isConvective = .FALSE.
                END IF
            END IF
            
            currentConv = intShell(ii)%isConvective
        END IF
        
        ! Check for lone radiative
        IF (currentConv.AND.(.NOT.nextConv)) THEN
            IF ((ii + 2).EQ.totShell) THEN
                intShell(ii + 1)%isConvective = .TRUE.
            ELSE IF (intShell(ii + 2)%isConvective) THEN
                intShell(ii + 1)%isConvective = .TRUE.
            END IF
        END IF
        
        ! Check for lone convective
        IF ((.NOT.currentConv).AND.nextConv) THEN
            IF ((ii + 2).EQ.totShell) THEN
                intShell(ii + 1)%isConvective = .FALSE.
            ELSE IF (.NOT.intShell(ii + 2)%isConvective) THEN
                intShell(ii + 1)%isConvective = .FALSE.
            END IF
        END IF
    END DO
    
END SUBROUTINE storeShells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine reduces the number of shells needed to integrate.        !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, inter-shell array.                                            !!!
!!! -liteShell, reduced inter-shell array.                                   !!!
!!! -extraShells, the c13 extra-shells array.                                !!!
!!! -totShell, size of model.                                                !!!
!!! -nLiteShell, size of liteShell.                                          !!!
!!! -siz, number of chemical species.                                        !!!
!!!                                                                          !!!
!!! On output, liteShell has the relevant intershell values.                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createLiteShells(intShell, liteShell, extraShells, totShell, &
                            nLiteShell, siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:), liteShell(:), extraShells(:)
    INTEGER::totShell, nLiteShell, siz
    
    ! Local
    DOUBLE PRECISION::extraMass0, extraMass1, massLimit = 5.D-7
    INTEGER::ii, lastII, extraSiz
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get extraShells limits
    IF (ASSOCIATED(extraShells)) THEN
        extraSiz = SIZE(extraShells)
        extraMass0 = extraShells(1)%mass0
        extraMass1 = extraShells(extraSiz)%mass1
    END IF
    
    ! Count how many shells should be in liteShell
    nLiteShell = 1; lastII = 1
    DO ii = 2, totShell - 1
        ! Count a shell only if greater than massLimit
        IF ((intShell(ii)%mass0 - intShell(lastII)%mass0).GE.massLimit) THEN
            ! Add count and advance lastII
            nLiteShell = nLiteShell + 1
            lastII = ii
        ELSE IF (ASSOCIATED(extraShells)) THEN
            ! Be careful with extraShells
            IF (intShell(ii)%mass0.GE.extraMass0) THEN
                IF (intShell(ii)%mass1.LE.extraMass1) THEN
                    ! Add count and advance lastII
                    nLiteShell = nLiteShell + 1
                    lastII = ii
                END IF
            END IF
        END IF
    END DO
    
    ! Allocate the liteShell
    ALLOCATE(liteShell(nLiteShell))
    DO ii = 1, nLiteShell
        ALLOCATE(liteShell(ii)%dens(siz))
    END DO
    
    ! Add first one
    liteShell(1) = intShell(1)
    liteShell(1)%dens = intShell(1)%dens
    
    ! Now fill them
    nLiteShell = 1; lastII = 1
    DO ii = 2, totShell - 1
        IF ((intShell(ii)%mass0 - intShell(lastII)%mass0).GE.massLimit) THEN
            ! Advance nLiteShell and add shell
            nLiteShell = nLiteShell + 1
            liteShell(nLiteShell) = intShell(ii)
            liteShell(nLiteShell)%dens = intShell(ii)%dens
            
            ! Advance lastII
            lastII = ii
        ELSE IF (ASSOCIATED(extraShells)) THEN
            ! Be careful with extraShells
            IF (intShell(ii)%mass0.GE.extraMass0) THEN
                IF (intShell(ii)%mass1.LE.extraMass1) THEN
                    ! Advance nLiteShell and add shell
                    nLiteShell = nLiteShell + 1
                    liteShell(nLiteShell) = intShell(ii)
                    liteShell(nLiteShell)%dens = intShell(ii)%dens
                    
                    ! Advance lastII
                    lastII = ii
                END IF
            END IF
        END IF
    END DO
    
END SUBROUTINE createLiteShells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine adds extra c13 shells.                                   !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, inter-shell array.                                            !!!
!!! -exShls, the extra shells array.                                         !!!
!!! -c13indx, corresponding index to c13.                                    !!!
!!! -p1indx, corresponding index to protons.                                 !!!
!!! -totShell, size of intShell + 1.                                         !!!
!!! -exNum, extra shells quantity.                                           !!!
!!! -lowMass, low mass limit on the extra shells.                            !!!
!!! -upMass, up mass limit on the extra shells.                              !!!
!!! -lastConv, true if last time we had a completely convective pocket.      !!!
!!! -yscale, threshold for zero abundances.                                  !!!
!!! -siz, number of chemical species.                                        !!!
!!!                                                                          !!!
!!! On output, intShell has exNum extra shells if the conditions are met.    !!!
!!! If conditions are not met, the exShls are cleaned if they existed.       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE addExtraShells(intShell, exShls, c13indx, p1indx, n14indx, &
                          totShell, exNum, lowMass, upMass, lastConv, yscale, &
                          siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:), exShls(:)
    DOUBLE PRECISION::lowMass, upMass, yscale
    INTEGER::c13indx, p1indx, n14indx, totShell, exNum, siz
    LOGICAL::lastConv
    
    ! Local
    TYPE (INTERSHELL), POINTER::storeIntShell(:)
    DOUBLE PRECISION, ALLOCATABLE::valDens(:, :), valMassDens(:), mass(:)
    DOUBLE PRECISION, ALLOCATABLE::interArr(:, :), interDens(:, :)
    DOUBLE PRECISION::valArray(14, totShell - 1), valMass(totShell - 1)
    DOUBLE PRECISION::dMass, holdMass, firstMass, lastMass, c13ShellVal, c13Lim
    DOUBLE PRECISION::p1Envelp, envFact, n14ShellVal
    INTEGER::ii, first, last, jj, kk, frstMin, lstMin, c13Mass, n14Mass
    INTEGER::cpExNum, preExNum, sizValMass
    LOGICAL::seenC13Pocket, doExtension
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Define c13Mass and n14Mass
    c13Mass = 13; n14Mass = 14
    
    ! Carbon limit in mass fraction
    c13Lim = 1.D-3
    
    ! Envelope factor (multiply evenlope hydrogen abundance by this
    ! to define threshold)
    envFact = 0.05
    
    ! Hydrogen envelope value
    p1Envelp = intShell(totShell - 1)%dens(p1indx)
    
    ! We want to enhance the c13 region that complies with 4 requierements:
    ! 1) It has to be greater than a certain value (c13Lim)
    ! 2) It has to be contiguous
    ! 3) It has to be in the intershell
    ! 4) Its appearance has to occur in the intershell-envelope transition
    
    ! We want to keep the extra points until the next pulse
    
    ! Check where c13 is contiguously greater than c13Lim
    seenC13Pocket = .FALSE.
    first = 0; last = 0; ii = 1
    IF (.NOT.ASSOCIATED(exShls)) THEN
        ! Look for pocket
        DO
            ! If in envelope, exit
            IF (intShell(ii)%dens(p1indx).GE.p1Envelp*envFact) EXIT
            
            ! Store value
            c13ShellVal = intShell(ii)%dens(c13indx)*c13Mass
            n14ShellVal = intShell(ii)%dens(n14indx)*n14Mass
            
            ! If seen C13, get pocket indices
            IF ((c13ShellVal.GT.c13Lim).AND.(c13ShellVal.GT.n14ShellVal)) THEN
                seenC13Pocket = .TRUE.
                first = ii
                DO
                    ! Store possible last value
                    last = ii
                    
                    ! Store value
                    c13ShellVal = intShell(ii)%dens(c13indx)*c13Mass
                    n14ShellVal = intShell(ii)%dens(n14indx)*n14Mass
                    
                    ! Check conditions
                    IF (c13ShellVal.LE.c13Lim) EXIT
                    IF (c13ShellVal.LE.n14ShellVal) EXIT
                    
                    ! Advance index
                    ii = ii + 1
                    IF (ii.GE.totShell) EXIT
                END DO
            END IF
            
            ! If we have the pocket, exit
            IF (seenC13Pocket) EXIT
            
            ! Advance index
            ii = ii + 1
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Store masses if identified pocket
        IF (seenC13Pocket) THEN
            lowMass = intShell(first)%mass0
            upMass = intShell(last)%mass1
            lastConv = .FALSE.
        END IF
    ELSE
        ! Search for indices and check if we can expand them
        
        ! Lower index
        DO ii = 1, totShell - 2
            IF (intShell(ii)%mass0.GT.lowMass) THEN
                first = ii - 1
                EXIT
            END IF
        END DO
        
        ! See if we can expand the lower index
        DO ii = first, 1, -1
            c13ShellVal = intShell(ii)%dens(c13indx)*c13Mass
            
            IF (c13ShellVal.LE.c13Lim) THEN
                ! Only change lowMass if found a possible expansion
                IF ((ii + 1).LT.first) THEN
                    first = ii + 1
                    lowMass = intShell(first)%mass0
                END IF
                
                EXIT
            END IF
        END DO
        
        ! Now upper index
        DO ii = totShell - 2, 1, -1
            IF (intShell(ii)%mass1.LT.upMass) THEN
                last = ii + 1
                EXIT
            END IF
        END DO
        
        ! See if we can expand the upper index
        DO ii = last, totShell - 2
            c13ShellVal = intShell(ii)%dens(c13indx)*c13Mass
            
            ! Don't expand if in envelope
            IF ((c13ShellVal.LE.c13Lim).OR. &
                (intShell(ii)%dens(p1indx).GT.p1Envelp*envFact)) THEN
                
                ! Only change upMass if found a possible expansion
                IF ((ii - 1).GT.last) THEN
                    last = ii - 1
                    upMass = intShell(last)%mass1
                END IF
                
                EXIT
            END IF
        END DO
    END IF
    
    ! Ensure dMass is not going to be too small
    IF (last.GT.0) THEN
        dMass = 1.D-14
        cpExNum = exNum + 1
        DO WHILE (dMass.LT.1.D-12)
            cpExNum = cpExNum - 1
            dMass = (intShell(last)%mass1 - intShell(first)%mass0)/cpExNum
        END DO
        
        ! If pocket has more points than extension, don't apply extension
        IF ((last - first).GT.cpExNum) seenC13Pocket = .FALSE.
    END IF
    
    ! Check if perform extension
    doExtension = (seenC13Pocket.OR.(ASSOCIATED(exShls).AND.(.NOT.(lastConv))))
    
    ! If from first to last everything is convective, then mark it
    IF (ASSOCIATED(exShls)) THEN
        lastConv = .TRUE.
        DO ii = first, last
            IF (.NOT.(intShell(ii)%isConvective)) THEN
                lastConv = .FALSE.
                EXIT
            END IF
        END DO
    END IF
    
    ! If everything is as we require, check if the exShls exist
    IF (doExtension) THEN
        IF (ASSOCIATED(exShls)) THEN
            ! Store previous exNum
            preExNum = SIZE(exShls)
            
            ! Check where do we need to introduce the exShls in intShell
            frstMin = 0; lstMin = 0
            firstMass = (exShls(1)%mass0 + exShls(1)%mass1)*0.5
            lastMass = (exShls(preExNum)%mass0 + exShls(preExNum)%mass1)*0.5
            DO ii = 1, totShell - 1
                IF (frstMin.EQ.0) THEN
                    ! Store masses
                    holdMass = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
                    
                    ! Compare them
                    IF (holdMass.GT.firstMass) THEN
                        ! Check that the previous one is lower, if it is not
                        ! we have to correct the index
                        holdMass = intShell(ii - 1)%mass0
                        holdMass = (holdMass + intShell(ii - 1)%mass1)*0.5
                        
                        ! Pick ii if previous mass is lower, otherwise (if it is
                        ! equal) pick ii - 1
                        IF (holdMass.LT.firstMass) THEN
                            frstMin = ii
                        ELSE
                            frstMin = ii - 1
                        END IF
                    END IF
                END IF
                
                IF (lstMin.EQ.0) THEN
                    ! Store masses
                    holdMass = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
                    
                    ! Compare them
                    IF (holdMass.GT.lastMass) THEN
                        lstMin = ii - 1
                        EXIT
                    END IF
                END IF
            END DO
            
            ! Allocate valMassDens and valDens
            ALLOCATE(valMassDens(totShell + preExNum - (lstMin - frstMin + 2)))
            ALLOCATE(valDens(siz, totShell + preExNum - (lstMin - frstMin + 2)))
            
            ! Fill them
            ii = 1; kk = 1
            sizValMass = SIZE(valMassDens)
            DO
                IF (ii.EQ.frstMin) THEN
                    DO jj = 1, preExNum
                        holdMass = (exShls(jj)%mass0 + exShls(jj)%mass1)*0.5
                        valMassDens(kk) = holdMass
                        valDens(:, kk) = exShls(jj)%dens
                        
                        kk = kk + 1
                    END DO
                    
                    ii = lstMin + 1
                END IF
                
                IF (kk.GT.sizValMass) EXIT
                
                holdMass = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
                valMassDens(kk) = holdMass
                valDens(:, kk) = intShell(ii)%dens
                
                kk = kk + 1; ii = ii + 1
                IF (ii.GE.totShell) EXIT
            END DO
            
            ! Deallocate exShls
            DO ii = 1, preExNum
                DEALLOCATE(exShls(ii)%dens)
            END DO
            DEALLOCATE(exShls); NULLIFY(exShls)
        ELSE
            ! Allocate valMassDens and valDens
            ALLOCATE(valMassDens(totShell - 1), valDens(siz, totShell - 1))
            
            ! Fill them
            DO ii = 1, totShell - 1
                valMassDens(ii) = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
                valDens(:, ii) = intShell(ii)%dens
            END DO
        END IF
        
        ! Create arrays for interpolation
        ! valMass, temp, rho, hp0, hp1, pres0, pres1, vel0, vel1, dens
        DO ii = 1, totShell - 1
            ! Masses
            valMass(ii) = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
            
            ! Physical values
            valArray(1, ii) = intShell(ii)%temp
            valArray(2, ii) = intShell(ii)%rho
            valArray(3, ii) = intShell(ii)%rho0
            valArray(4, ii) = intShell(ii)%rho1
            valArray(5, ii) = intShell(ii)%hp0
            valArray(6, ii) = intShell(ii)%hp1
            valArray(7, ii) = intShell(ii)%pres0
            valArray(8, ii) = intShell(ii)%pres1
            valArray(9, ii) = intShell(ii)%vel0
            valArray(10, ii) = intShell(ii)%vel1
            valArray(11, ii) = intShell(ii)%rad0
            valArray(12, ii) = intShell(ii)%rad1
            valArray(13, ii) = intShell(ii)%radiat0
            valArray(14, ii) = intShell(ii)%radiat1
        END DO
        
        ! Create extra shells
        ALLOCATE(exShls(cpExNum))
        DO ii = 1, cpExNum
            ALLOCATE(exShls(ii)%dens(siz))
        END DO
        
        ! Allocate needed arrays
        ALLOCATE(interArr(SIZE(valArray(:, 1)), cpExNum))
        ALLOCATE(interDens(siz, cpExNum), mass(cpExNum))
        
        ! Mass array
        dMass = (intShell(last)%mass1 - intShell(first)%mass0)/cpExNum
        mass(1) = intShell(first)%mass0 + dMass*0.5
        DO ii = 2, cpExNum
            mass(ii) = mass(ii - 1) + dMass
        END DO
        
        ! Interpolate
        DO ii = 1, SIZE(interArr(:, 1))
            interArr(ii, :) = interpolate(mass, valMass, valArray(ii, :), &
                                          cpExNum)
        END DO
        DO ii = 1, siz
            interDens(ii, :) = interpolate(mass, valMassDens, valDens(ii, :), &
                                           cpExNum)
        END DO
        
        ! Put values in the extra shells
        exShls(1)%mass0 = intShell(first)%mass0
        exShls(1)%mass1 = exShls(1)%mass0 + dMass
        DO ii = 1, cpExNum
            IF (ii.NE.1) THEN
                exShls(ii)%mass0 = exShls(ii - 1)%mass1
                exShls(ii)%mass1 = exShls(ii)%mass0 + dMass
            END IF
            
            exShls(ii)%temp = interArr(1, ii)
            exShls(ii)%rho = interArr(2, ii)
            exShls(ii)%rho0 = interArr(3, ii)
            exShls(ii)%rho1 = interArr(4, ii)
            exShls(ii)%hp0 = interArr(5, ii)
            exShls(ii)%hp1 = interArr(6, ii)
            exShls(ii)%pres0 = interArr(7, ii)
            exShls(ii)%pres1 = interArr(8, ii)
            exShls(ii)%vel0 = interArr(9, ii)
            exShls(ii)%vel1 = interArr(10, ii)
            exShls(ii)%rad0 = interArr(11, ii)
            exShls(ii)%rad1 = interArr(12, ii)
            exShls(ii)%radiat0 = interArr(13, ii)
            exShls(ii)%radiat1 = interArr(14, ii)
            
            exShls(ii)%dens = interDens(:, ii)
        END DO
        
        ! Fix their convective status
        jj = 1
        DO ii = first, MIN(last + 1, totShell - 1)
            DO WHILE (exShls(jj)%mass1.LE.intShell(ii)%mass1)
                exShls(jj)%isConvective = intShell(ii)%isConvective
                exShls(jj)%convShell => intShell(ii)%convShell
                
                jj = jj + 1
                IF (jj.GT.cpExNum) EXIT
            END DO
            
            IF (jj.GT.cpExNum) EXIT
        END DO
        
        ! Add the extra shells to intShell
        
        ! Copy and nullify
        storeIntShell => intShell
        NULLIFY(intShell)
        
        ! Allocate
        ALLOCATE(intShell(totShell + cpExNum - (last - first + 2)))
        
        ! Add them
        intShell(1:first - 1) = storeIntShell(1:first - 1)
        intShell(first:first + cpExNum - 1) = exShls
        intShell(first + cpExNum:) = storeIntShell(last + 1:)
        
        ! Free memory
        DEALLOCATE(interArr, interDens, mass, valMassDens, valDens)
        
        DO ii = 1, SIZE(storeIntShell)
            DEALLOCATE(storeIntShell(ii)%dens)
        END DO
        DEALLOCATE(storeIntShell); NULLIFY(storeIntShell)
        
        ! Update totShell
        totShell = SIZE(intShell) + 1
    ELSE
        ! Deallocate the extra shells if allocated
        IF (ASSOCIATED(exShls)) THEN
            DO ii = 1, cpExNum
                DEALLOCATE(exShls(ii)%dens)
            END DO
            DEALLOCATE(exShls); NULLIFY(exShls)
        END IF
    END IF
    
    IF (doExtension) THEN
        ! Check if negatives exist and make them 0
        DO ii = 1, totShell - 1
            DO kk = 1, siz
                IF (ABS(intShell(ii)%dens(kk)).LT.yscale) THEN
                    intShell(ii)%dens(kk) = 0.D0
                ELSE IF (intShell(ii)%dens(kk).LT.0.D0) THEN
                    PRINT*, "# Extra shells routine"
                    PRINT*, "# Warning, making zero shell and element ", ii, kk
                    PRINT*, "# Previous value: ", intShell(ii)%dens(kk)
                    intShell(ii)%dens(kk) = 0.D0
                END IF
            END DO
        END DO
    END IF
    
END SUBROUTINE addExtraShells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine interpolates liteShell into intShell.                    !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, inter-shell array.                                            !!!
!!! -liteShell, reduced inter-shell array.                                   !!!
!!! -siz, number of chemical species.                                        !!!
!!!                                                                          !!!
!!! On output, intShell has the interpolated values from liteShell.          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE undoLiteShell(intShell, liteShell, siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:), liteShell(:)
    INTEGER::siz
    
    ! Local
    DOUBLE PRECISION, POINTER::mass(:), valMass(:), valDens(:, :), dens(:)
    INTEGER::ii, jj, nIntShell, nLiteShell
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get values
    nIntShell = SIZE(intShell)
    nLiteShell = SIZE(liteShell)
    ALLOCATE(dens(nIntShell))
    ALLOCATE(mass(nIntShell), valMass(nLiteShell), valDens(nLiteShell, siz))
    
    ! Create the mass, valMass and valDens arrays
    DO ii = 1, nIntShell
        mass(ii) = (intShell(ii)%mass0 + intShell(ii)%mass1)*0.5
    END DO
    DO ii = 1, nLiteShell
        valMass(ii) = (liteShell(ii)%mass0 + liteShell(ii)%mass1)*0.5
        valDens(ii, :) = liteShell(ii)%dens
    END DO
    
    ! Interpolate each species
    DO ii = 1, siz
        dens = ABS(interpolate(mass, valMass, valDens(:, ii), nIntShell))
        DO jj = 1, nIntShell
            intShell(jj)%dens(ii) = dens(jj)
        END DO
    END DO
    
    ! Deallocate pointers
    DEALLOCATE(mass, valMass, valDens, dens)
    
END SUBROUTINE undoLiteShell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine undoes the shell transformations.                        !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, inter-shell array.                                            !!!
!!! -model, array which will hold the updated chemistry.                     !!!
!!! -totShell, size of model.                                                !!!
!!!                                                                          !!!
!!! On output, model has the correct values from intShell.                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE undoTransformations(intShell, model, totShell)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:)
    TYPE (SHELL), POINTER::model(:)
    INTEGER::totShell
    
    ! Local
    TYPE (SHELL), POINTER::tempModel(:)
    INTEGER::ii, jj, sizModel
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Copy original model
    sizModel = SIZE(model)
    ALLOCATE(tempModel(sizModel))
    DO ii = 1, sizModel
        tempModel(ii)%temp = model(ii)%temp
        tempModel(ii)%rho = model(ii)%rho
        tempModel(ii)%mass = model(ii)%mass
        tempModel(ii)%radius = model(ii)%radius
        tempModel(ii)%hp = model(ii)%hp
        tempModel(ii)%pressure = model(ii)%pressure
        tempModel(ii)%velocity = model(ii)%velocity
        tempModel(ii)%radiat = model(ii)%radiat
    END DO
    
    ! Deallocate and reallocate the model
    DO ii = 1, sizModel
        DEALLOCATE(model(ii)%dens)
    END DO
    DEALLOCATE(model)
    
    ALLOCATE(model(totShell))
    
    ! Introduce every "dens" value and only the new physical values
    model(1)%dens = intShell(1)%dens
    model(1)%temp = intShell(1)%temp
    model(1)%rho = intShell(1)%rho
    
    jj = 1 ! Index for tempModel
    DO ii = 1, totShell - 1
        ! Abundance, temp and rho
        ! Temp and rho have to be written always in intershell form
        IF (ii.GT.1) THEN
            model(ii)%dens = 2*intShell(ii - 1)%dens - model(ii - 1)%dens
            model(ii)%temp = 2*intShell(ii - 1)%temp - model(ii - 1)%temp
            model(ii)%rho = 2*intShell(ii - 1)%rho - model(ii - 1)%rho
        END IF
        
        ! Advance index if bigger mass in intShell
        DO WHILE ((intShell(ii)%mass0 - tempModel(jj)%mass).GT.1.D-14)
            jj = jj + 1
        END DO
        
        ! Physics
        IF (ABS(intShell(ii)%mass0 - tempModel(jj)%mass).LT.1.D-14) THEN
            ! Pure shell values
            model(ii)%mass = tempModel(jj)%mass
            model(ii)%radius = tempModel(jj)%radius
            model(ii)%hp = tempModel(jj)%hp
            model(ii)%pressure = tempModel(jj)%pressure
            model(ii)%velocity = tempModel(jj)%velocity
            model(ii)%radiat = tempModel(jj)%radiat
            
            jj = jj + 1
        ELSE
            ! Safeguard
            IF (ii.EQ.totShell - 1) THEN
                PRINT*, "WARNING: In undoTransformations"
                PRINT*, "ii.EQ.totShell - 1 for intershell values"
                STOP
            END IF
            
            ! Intershell values
            model(ii)%mass = intShell(ii)%mass0
            model(ii)%radius = intShell(ii)%rad0
            model(ii)%hp = intShell(ii)%hp0
            model(ii)%pressure = intShell(ii)%pres0
            model(ii)%velocity = intShell(ii)%vel0
            model(ii)%radiat = intShell(ii)%radiat0
        END IF
    END DO
    
    ! Suppose last values are pure shell values
    model(totShell)%dens = 2*intShell(totShell - 1)%dens
    model(totShell)%dens = model(totShell)%dens - model(totShell - 1)%dens
    
    model(totShell)%temp = 2*intShell(totShell - 1)%temp
    model(totShell)%temp = model(totShell)%temp - model(totShell - 1)%temp
    
    model(totShell)%rho = 2*intShell(totShell - 1)%rho
    model(totShell)%rho = model(totShell)%rho - model(totShell - 1)%rho
    
    model(totShell)%mass = tempModel(sizModel)%mass
    model(totShell)%radius = tempModel(sizModel)%radius
    model(totShell)%hp = tempModel(sizModel)%hp
    model(totShell)%pressure = tempModel(sizModel)%pressure
    model(totShell)%velocity = tempModel(sizModel)%velocity
    model(totShell)%radiat = tempModel(sizModel)%radiat
    
    ! Deallocate tempModel
    DEALLOCATE(tempModel)
    
END SUBROUTINE undoTransformations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! The next subroutine clears the reactions linked list.                    !!!
!!!                                                                          !!!
!!! The input value is:                                                      !!!
!!! -node, a REACT type pointer.                                             !!!
!!!                                                                          !!!
!!! On output, the list has been cleared from tail to head.                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECURSIVE SUBROUTINE clearNodes(node)
    IMPLICIT NONE
    
    ! Input
    TYPE (REACT), POINTER::node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF (ASSOCIATED(node%next)) CALL clearNodes(node%next)
    IF (ALLOCATED(node%crossTable)) DEALLOCATE(node%crossTable)
    IF (ALLOCATED(node%avector)) DEALLOCATE(node%avector)
    
    DEALLOCATE(node)
    NULLIFY(node)
    
END SUBROUTINE clearNodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! The next subroutine clears the cross sections linked list.               !!!
!!!                                                                          !!!
!!! The input value is:                                                      !!!
!!! -node, a CROSSARR type pointer.                                          !!!
!!!                                                                          !!!
!!! On output, the list has been cleared from tail to head.                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECURSIVE SUBROUTINE clearCrossNodes(node)
    IMPLICIT NONE
    
    ! Input
    TYPE (CROSSARR), POINTER::node

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IF (ASSOCIATED(node%next)) CALL clearCrossNodes(node%next)
    
    DEALLOCATE(node)
    NULLIFY(node)
    
END SUBROUTINE clearCrossNodes

END MODULE shell_operations

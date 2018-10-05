MODULE integrator
    USE MPI
    USE integvars_mod
    USE structures_mod
    USE math_routines
    USE mixer
    IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine performs a coupled chemical and mixing step.             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -dt, the total timestep for integrating in seconds.                      !!!
!!! -intShell, the intershells.                                              !!!
!!! -totShell, number of intershells +1.                                     !!!
!!! -crosLst, the array with all cross sections stored.                      !!!
!!! -ntwkMass, an array with the atomic weights of the species in order.     !!!
!!! -eIndices, an array of type INDICES with the indices of the maximum      !!!
!!!             contributors to the free electrons in the star.              !!!
!!! -p1indx, index at which protons are.                                     !!!
!!! -he4indx, index at which he4 is.                                         !!!
!!! -ovParam, the envelope overshooting parameter.                           !!!
!!! -ovPDCZParam, the PDCZ overshooting parameter.                           !!!
!!! -mixFreq, frequency of mixing in the global timestep.                    !!!
!!! -eps, the relative accuracy.                                             !!!
!!! -siz, the number of species (size of dens).                              !!!
!!! -yscale, the threshold below which the abundances are considered zero.   !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -firstIntegShell, first shell to consider integration.                   !!!
!!! -lastIntegShell, last shell to consider integration.                     !!!
!!! -ovMode, the overshooting mode (advective or diffusive).                 !!!
!!! -nProc, total number of threads.                                         !!!
!!! -rank, this thread index.                                                !!!
!!!                                                                          !!!
!!! At the output, intShell%dens has the updated abundances.                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mixedIntegration(dt, intShell, totShell, crosLst, ntwkMass, &
                            eIndices, p1indx, he4indx, ovParam, ovPDCZParam, &
                            mixFreq, eps, siz, yscale, firstOv, &
                            firstIntegShell, lastIntegShell, ovMode, nProc, &
                            rank)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL)::intShell(:)
    TYPE (CROSSARR)::crosLst(:)
    TYPE (INDICES)::eIndices(:)
    DOUBLE PRECISION::dt, ovParam, ovPDCZParam, eps, yscale
    INTEGER::totShell, ntwkMass(:), p1indx, he4indx, mixFreq, siz
    INTEGER::firstIntegShell, lastIntegShell, firstOv, nProc, rank
    CHARACTER(20)::ovMode
    
    ! Local
    TYPE (INTERSHELL), POINTER::ovShell(:)
    DOUBLE PRECISION, POINTER::ovMatrix(:, :), dx(:)
    DOUBLE PRECISION::mixHH
    INTEGER, POINTER::convecIndex(:, :)
    INTEGER::firstMix, ii, kk, nShells
    LOGICAL::performOv
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Check if negative values
    DO ii = 1, totShell - 1
        IF (MINVAL(intShell(ii)%dens).LT.0.D0) THEN
            PRINT*, "Initial negative value!"
            PRINT*, "Indices:"
            
            DO kk = 1, siz
                IF (intShell(ii)%dens(kk).LT.0.D0) THEN
                    PRINT*, kk, intShell(ii)%dens(kk)
                END IF
            END DO
            
            STOP
        END IF
    END DO
    
    ! Mix convection before starting
    CALL mixConvection(intShell, totShell, siz)
    
    ! Create Overshooting Matrix
    IF ((ovParam.GT.0.D0).OR.(ovPDCZParam.GT.0.D0)) THEN
        IF (ovMode.EQ."advective") THEN
            CALL createOvMatrix(ovMatrix, intShell, ovShell, nShells, firstOv, &
                                p1indx, he4indx, ovParam, ovPDCZParam, &
                                convecIndex, performOv)
        ELSE IF (ovMode.EQ."diffusive") THEN
            CALL createOvArrays(ovMatrix, intShell, ovShell, nShells, firstOv, &
                                p1indx, he4indx, ovParam, ovPDCZParam, &
                                convecIndex, dx, performOv)
        END IF
    ELSE
        performOv = .FALSE.
    END IF
    
    ! Get first mixed shell
    firstMix = lastIntegShell
    
    ! Overshooting
    IF (performOv) firstMix = firstOv
    
    ! Convection
    DO ii = firstIntegShell, lastIntegShell
        IF (intShell(ii)%isConvective) THEN
            IF (ii.LT.firstMix) firstMix = ii
            EXIT
        END IF
    END DO
    
    ! Integrate non-mixed shells
    CALL noMixIntegration(dt, intShell, totShell, crosLst, ntwkMass, eIndices, &
                        eps, siz, yscale, firstIntegShell, firstMix - 1, &
                        nProc, rank)
    
    ! Now integrate and mix as many times as needed by mixFreq
    mixHH = dt/mixFreq
    DO ii = 1, mixFreq
        CALL noMixIntegration(mixHH, intShell, totShell, crosLst, ntwkMass, &
                            eIndices, eps, siz, yscale, firstMix, &
                            lastIntegShell, nProc, rank)
        
        ! Now mix
        CALL mixConvection(intShell, totShell, siz)
        
        ! Check if performing overshooting
        IF (.NOT.performOv) CYCLE
        
        CALL applyOvershooting(intShell, totShell, ovMatrix, ovShell, nShells, &
                        firstOv, mixHH, eps, convecIndex, yscale, siz, dx, &
                        ovMode, nProc, rank)
        
        CALL mixConvection(intShell, totShell, siz)
    END DO
    
    ! Deallocate ovMatrix, ovShell and "dens" in shellCpy
    IF ((ovParam.GT.0.D0).OR.(ovPDCZParam.GT.0.D0)) THEN
        DO ii = 1, SIZE(ovShell)
            DEALLOCATE(ovShell(ii)%dens)
        END DO
        DEALLOCATE(ovMatrix, ovShell)
    END IF
    
    ! Deallocate convecIndex
    IF (performOv.OR.(ovMode.EQ."diffusive")) DEALLOCATE(convecIndex)
    IF (ovMode.EQ."diffusive") DEALLOCATE(dx)
    
END SUBROUTINE mixedIntegration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine integrates without mixing.                               !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -dt, the total timestep for integrating in seconds.                      !!!
!!! -intShell, the intershells.                                              !!!
!!! -totShell, number of intershells +1.                                     !!!
!!! -crosLst, the array with all cross sections stored.                      !!!
!!! -ntwkMass, an array with the atomic weights of the species in order.     !!!
!!! -eIndices, an array of type INDICES with the indices of the maximum      !!!
!!!             contributors to the free electrons in the star.              !!!
!!! -eps, the relative accuracy.                                             !!!
!!! -siz, the number of species (size of dens).                              !!!
!!! -yscale, the threshold below which the abundances are considered zero.   !!!
!!! -firstIntegShell, first shell to consider integration.                   !!!
!!! -lastIntegShell, last shell to consider integration.                     !!!
!!! -nProc, total number of threads.                                         !!!
!!! -rank, this thread index.                                                !!!
!!!                                                                          !!!
!!! At the output, intShell%dens has the updated abundances.                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE noMixIntegration(dt, intShell, totShell, crosLst, ntwkMass, &
                        eIndices, eps, siz, yscale, firstIntegShell, &
                        lastIntegShell, nProc, rank)
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::ierror
    
    ! Input
    TYPE (INTERSHELL)::intShell(:)
    TYPE (CROSSARR)::crosLst(:)
    TYPE (INDICES)::eIndices(:)
    DOUBLE PRECISION::dt, eps, yscale
    INTEGER::totShell, ntwkMass(:), siz, nProc, rank
    INTEGER::firstIntegShell, lastIntegShell
    
    ! Local
    DOUBLE PRECISION::dens(siz, totShell - 1), initialMass(totShell - 1), Htot
    DOUBLE PRECISION::tt(siz, nSubsteps, nSubsteps), htm, err, epsFactor = 0.95
    DOUBLE PRECISION::error(nSubsteps - 1), ym(siz), hhcoef, delt(siz)
    DOUBLE PRECISION::newHH, difference, minHH, HH, work(nSubsteps - 1), minWork
    INTEGER::ii, mm, kk, lastm, mink, calcProc
    LOGICAL::skipStep, converged, firstTry
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Use the explicit method for every shell
    
    ! Copy densities, initialize initial mass
    DO ii = 1, totShell - 1
        dens(:, ii) = intShell(ii)%dens
        initialMass(ii) = SUM(dens(:, ii)*ntwkMass)
    END DO
    
    ! Initialize variables outside the loop
    minHH = dt
    
    ! Define aj
    IF (aj(1).EQ.1) THEN
        DO ii = 1, SIZE(substeps)
            aj(ii) = SUM(substeps(1:ii))
        END DO
        
        DO ii = 1, SIZE(substeps) - 1
            DO kk = 1, ii - 1
                alph(kk, ii) = (aj(kk + 1) - aj(ii + 1))
                alph(kk, ii) = alph(kk, ii)/((kk + 1)*(aj(ii + 1) - aj(1) + 1))
                alph(kk, ii) = eps**alph(kk, ii)
            END DO
        END DO
    END IF
    
    ! Main loop
    DO ii = firstIntegShell, lastIntegShell
        IF (nProc.GT.1) THEN
            calcProc = MOD(ii - 1, nProc)
            IF (calcProc.NE.rank) CYCLE
        END IF
        
        ! Initialize variables
        converged = .TRUE.; firstTry = .TRUE.
        Htot = dt; HH = minHH; newHH = HH
        
        DO WHILE(Htot.GT.0.D0)
            skipStep = .FALSE.
            tt = 0.D0; error = 0.D0; lastm = 1
            
            ! Copy new abundances
            intShell(ii)%dens = dens(:, ii)
            
            ! This loop performs the actual integration with the HH step.
            DO mm = 1, nSubsteps
                lastm = mm
                htm = HH/substeps(mm)
                
                ! Calculate the deltas and "final" densities for this step.
                ym = intShell(ii)%dens; delt = 0.D0
                DO kk = 1, substeps(mm)
                    CALL solveExplicit(htm, ym, crosLst(ii), siz, eIndices, &
                                       skipStep)
                    
                    IF (skipStep) EXIT
                END DO
                
                ! Exit if skipped
                IF (skipStep) EXIT
                
                ! Extrapolate.
                tt(:, 1, mm) = ym
                DO kk = 2, mm
                    tt(:, kk, mm) = tt(:, kk - 1, mm) + &
                            (tt(:, kk - 1, mm) - tt(:, kk - 1, mm - 1))/&
                            ((FLOAT(substeps(mm))/substeps(mm - kk + 1)) - 1)
                END DO
                
                ! Calculate error.
                IF (mm.GT.1) THEN
                    error(mm - 1) = 0.D0
                    
                    DO kk = 1, siz
                        ! Add this error
                        err = tt(kk, mm - 1, mm) - tt(kk, mm, mm)
                        
                        IF (ABS(tt(kk, mm - 1, mm)).GE.yscale) THEN
                            err = (err/tt(kk, mm - 1, mm))**2
                        ELSE
                            err = (err/yscale)**2
                        END IF
                        
                        error(mm - 1) = error(mm - 1) + err
                    END DO
                    
                    ! Calculate the average
                    error(mm - 1) = SQRT(error(mm - 1)/siz)
                    
                    ! If error < eps, we don't have to continue calculating
                    IF (error(mm - 1).LT.eps) EXIT
                END IF
            END DO
            
            ! If skip step calculate HH and return to the beginning
            IF (skipStep) THEN
                ! Skipped step, reduce fast
                IF (ABS(HH - 1.D0).LT.0.5) THEN
                    HH = HH*0.5
                ELSE IF (HH.GT.1.D0) THEN
                    HH = SQRT(HH)
                ELSE
                    HH = HH**2
                END IF
                
                CYCLE
            END IF
            
            ! Apply absolute value to solution, decreasing the error:
            ! Exact solution A, approximate solution A' = A + p(h) where
            ! h is the step and p is a polynomial. We know A >= 0 and in that
            ! case ||A + p(h)| - A| <= |p(h)|
            !
            ! In fact, in general if A' < 0 and A > 0, if you take A'' = -k*A'
            ! then |A'' - A| < |A' - A| if -1 < k < 1 - 2*A/A' with the special
            ! case of |A'' - A| = 0 if k = -A/A'
            tt(:, lastm - 1, lastm) = ABS(tt(:, lastm - 1, lastm))
            
            ! Check error.
            converged = .FALSE.
            difference = SUM(tt(:, lastm - 1, lastm)*ntwkMass) - initialMass(ii)
            IF (ABS(difference).GT.(eps*epsFactor)) THEN
                ! If the difference is lower than eps, then update initialMass
                IF (ABS(difference).LE.eps) THEN
                    initialMass(ii) = SUM(dens(:, ii)*ntwkMass)
                END IF
            ELSE IF (error(lastm - 1).LT.eps) THEN
                intShell(ii)%dens = tt(:, lastm - 1, lastm)
                converged = .TRUE.
            END IF
            
            ! Advance solution if converged
            IF (converged) THEN
                Htot = Htot - HH
                
                ! Update abundances and initial mass for next step
                dens(:, ii) = intShell(ii)%dens
                initialMass(ii) = SUM(dens(:, ii)*ntwkMass)
            END IF
            
            ! New HH:
            ! First the root error and work.
            DO kk = 1, lastm - 1
                error(kk) = (error(kk)/(epsFactor*eps))**(1.D0/(kk + 2))
                work(kk) = error(kk)*aj(kk + 1)
            END DO
            
            ! Pick the minimum work
            minWork = work(1); mink = 1
            DO kk = 2, lastm - 1
                IF (work(kk).LT.minWork) THEN
                    minWork = work(kk)
                    mink = kk
                END IF
            END DO
            
            ! Reduce only if did not converge
            IF (.NOT.converged) THEN
                hhcoef = error(mink)
            ELSE
                hhcoef = error(lastm - 1)
                mink = lastm - 1
            END IF
            
            ! Increase order
            IF (((mink + 1).EQ.lastm).AND.((mink + 1).LT.nSubsteps)) THEN
                hhcoef = hhcoef/alph(lastm, lastm + 1)
            ELSE
                DO kk = 1, mink
                    IF ((hhcoef*alph(kk, mink + 2)).LT.error(kk)) THEN
                        hhcoef = error(kk)/alph(kk, mink + 1)
                        EXIT
                    END IF
                END DO
            END IF
            
            ! Now calculate the new HH (speeding things up a bit):
            IF (hhcoef.LE.0.D0) THEN
                newHH = Htot
            ELSE IF (.NOT.converged) THEN
                newHH = HH*1.D-1
            ELSE
                ! Adjust minHH
                IF (HH.LT.minHH) minHH = newHH
                IF (firstTry) THEN
                    IF (lastm.LE.5) minHH = dt
                    firstTry = .FALSE.
                END IF
                
                newHH = HH/hhcoef
            END IF
            
            HH = MINVAL((/newHH, Htot/))
            
            ! If HH is a 10% lower or less than the remaining Htot,
            ! take Htot instead
            IF ((HH*1.1).GE.Htot) HH = Htot
        END DO
    END DO
    
    ! Broadcast "dens" values:
    IF (nProc.GT.1) THEN
        DO ii = 1, totShell - 1
            ! Calculate processor rank
            calcProc = MOD(ii - 1, nProc)
            
            ! Broadcast
            CALL MPI_BCAST(dens(:, ii), siz, MPI_DOUBLE_PRECISION, calcProc, &
                           MPI_COMM_WORLD, ierror)
        END DO
    END IF
    
    ! Update final abundances
    DO ii = 1, totShell - 1
        intShell(ii)%dens = dens(:, ii)
    END DO
    
END SUBROUTINE noMixIntegration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine advances the calculation with the explicit method.       !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -htm, the current precision order timestep.                              !!!
!!! -dens, the abundances in molar fraction [mol/g].                         !!!
!!! -crosLst, a linked list of cross sections correspondent with reacts.     !!!
!!! -siz, the number of species (size of dens).                              !!!
!!! -eIndices, an array of type INDICES with the indices of the maximum      !!!
!!!            contributors to the free electrons in the star.               !!!
!!! -skipStep, boolean value which returns .T. if numerical divergence.      !!!
!!!                                                                          !!!
!!! At the output, dens is modified with the partial integration and         !!!
!!! skipStep might be .T. if there was a divergence.                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solveExplicit(htm, dens, crosLst, siz, eIndices, skipStep)
    IMPLICIT NONE
    
    ! Input
    TYPE (CROSSARR)::crosLst
    TYPE (INDICES)::eIndices(:)
    DOUBLE PRECISION::htm, dens(:)
    INTEGER::siz
    LOGICAL::skipStep
    
    ! Local
    DOUBLE PRECISION::kk(siz), dd(siz)
    INTEGER::ii
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Calculate the rates and integrate
    CALL createRates(crosLst, dens, kk, dd, eIndices)
    dens = (dens + kk*htm)/(1 + dd*htm)
    
    ! Check current ym for big values or NaN
    DO ii = 1, siz
        IF ((ABS(dens(ii)).GT.1.D10).OR.(ISNAN(dens(ii)))) THEN
            skipStep = .TRUE.
            EXIT
        END IF
    END DO
    
END SUBROUTINE solveExplicit

END MODULE integrator

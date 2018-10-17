!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                                                                          !!!
!!! Slow process NUcleosynthesis Post-Processing code for ATon: SNUPPAT      !!!
!!!                                                                          !!!
!!! Written by Andrés Yagüe López in collaboration with Dr. Paolo Ventura,   !!!
!!! Dr. Aníbal García Hernández and Dr. Maria Lugaro.                        !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM snuppat
    USE MPI
    USE structures_mod
    USE loader
    USE integrator
    USE writer_reader
    USE mixer
    USE math_routines
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::rank, nProc, ierror
    
    ! Program variables
    TYPE (REACT), POINTER::reacts, lowTempReacts, highTempReacts
    TYPE (CROSSARR), POINTER::shellReacts(:)
    TYPE (SHELL), POINTER::model(:), nextModel(:)
    TYPE (INTERSHELL), POINTER::intShell(:), extraShells(:), liteShell(:)
    TYPE (INDICES)::eIndices(5), capIndices(7)
    DOUBLE PRECISION, ALLOCATABLE::partfunct(:, :), dens(:)
    DOUBLE PRECISION::t1, t2, mass, time1, time2, eps, temp, rho, tTot, modDt
    DOUBLE PRECISION::yscale, ovParam, sInYear = 31556926.D0, intDt, intT1
    DOUBLE PRECISION::minIntDt, minPulsIntDt, lowMass, upMass, ovPDCZParam
    DOUBLE PRECISION::metallicity
    INTEGER, ALLOCATABLE::ntwkMass(:)
    INTEGER::siz, totShell, ii, jj, jjindx, intgShls, modNum, newModNum, n1indx
    INTEGER::writeFreq, nSteps, p1indx, he4indx, c13indx, n14indx, mixFreq
    INTEGER::firstIntegShell, lastIntegShell, calcProc, extraNum, firstOv
    INTEGER::nLiteShell
    CHARACTER(20)::ovMode
    LOGICAL::integThis, isNext, isLowTemp, isMaster, isPulse, performDt, cont
    LOGICAL::lastConv
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize MPI
    CALL MPI_INIT(ierror)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nProc, ierror)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
    
    ! Identify master
    isMaster = (rank.EQ.0)
    
    ! Read simulation parameters
    OPEN(UNIT = uni, FILE = "parameters.in")
    READ(uni, *) eps
    READ(uni, *) yscale
    READ(uni, *) extraNum
    READ(uni, *) ovParam
    READ(uni, *) ovPDCZParam
    READ(uni, *) ovMode
    READ(uni, *) metallicity
    READ(uni, *) mixFreq
    READ(uni, *) minIntDt
    READ(uni, *) minPulsIntDt
    READ(uni, *) writeFreq
    CLOSE(UNIT = uni)
    
    ! Check that ovMode has one of the two correct values
    IF ((ovMode.NE."advective").AND.(ovMode.NE."diffusive")) THEN
        IF (isMaster) THEN
            PRINT*, "ovMode should be either 'advective' or 'diffusive'"
        END IF
        CALL MPI_FINALIZE(ierror)
        STOP
    END IF
    
    ! Open list with extra information
    OPEN(UNIT = uni, FILE = "data/"//filenames(1)//".lst")
    
    ! Skip first line
    READ(uni, *)
    
    ! Store in siz the number of elements.
    READ(uni, *) siz
    
    ! Read eIndices
    DO ii = 1, SIZE(eIndices)
        READ(uni, *) eIndices(ii)%elem
        READ(uni, *) eIndices(ii)%posit
    END DO
    
    ! Read capIndices
    ! Here we have in this order: neutrons, protons, he4, c13, n14, ne22
    DO ii = 1, SIZE(capIndices)
        READ(uni, *) capIndices(ii)%elem
        READ(uni, *) capIndices(ii)%posit
    END DO
    
    CLOSE(UNIT = uni)
    
    ! Store cap indices positions
    DO ii = 1, SIZE(capIndices)
        IF (capIndices(ii)%elem.EQ."n") n1indx = capIndices(ii)%posit
        IF (capIndices(ii)%elem.EQ."p") p1indx = capIndices(ii)%posit
        IF (capIndices(ii)%elem.EQ."he4") he4indx = capIndices(ii)%posit
        IF (capIndices(ii)%elem.EQ."c13") c13indx = capIndices(ii)%posit
        IF (capIndices(ii)%elem.EQ."n14") n14indx = capIndices(ii)%posit
    END DO
    
    ! Allocate known arrays
    ALLOCATE(partfunct(siz, 24), dens(siz), ntwkMass(siz))
    
    ! Nullify non-allocated arrays
    NULLIFY(extraShells, reacts, model, nextModel, intShell)
    NULLIFY(lowTempReacts, highTempReacts)
    
    ! Signal start
    IF (isMaster) THEN
        PRINT*, "# Starting simulation"
        PRINT*, "# Overshooting parameters: ", ovParam, ovPDCZParam
        
        ! Load network and partition functions as they are immutable
        PRINT*, "# Loading network"
    END IF
    
    partfunct = 0.D0
    CALL loadNetwork(highTempReacts, lowTempReacts)
    CALL loadPartitions(partfunct)
    
    ! Open physical file
    OPEN(UNIT = uni, FILE = physics)
    
    ! Prepare first model
    IF (isMaster) PRINT*, "# Preparing first model"
    cont = .FALSE.
    CALL initialModel(model, t1, ntwkMass, modNum, metallicity, n1indx, &
                      isPulse, siz, cont, rank)
    totShell = SIZE(model)
    
    ! Output blank division for actual integration
    IF (isMaster) PRINT*
    
    ! Initialize
    nSteps = 0; tTot = 0.D0; intT1 = t1
    
    ! Make sure to always write the first model
    IF (.NOT.cont) nSteps = writeFreq - 1
    
    ! Integrate
    performDt = .TRUE.
    DO
        CALL CPU_TIME(time1)
        
        ! Read next model (we need dt)
        isNext = readPhysMod(nextModel, mass, t2, newModNum, isPulse, 0, siz)
        
        IF (.NOT.isNext) THEN
            IF (isMaster) THEN
                PRINT*, "# End of physical models"
                PRINT*, "# Overshooting parameter used: ", ovParam, ovPDCZParam
            END IF
            
            EXIT
        END IF
        
        ! This model dt in years
        modDt = (1.D1**t2 - 1.D1**t1)
        
        IF (isMaster) THEN
            PRINT*, "# Model number: ", newModNum
            PRINT*, "# Timestep years: ", modDt
        END IF
        
        ! Perform overshooting and integration
        ! ====================================
        intgShls = 0
        
        IF (performDt) THEN
            ! Allocate intShell
            ALLOCATE(intShell(totShell - 1))
            
            ! Store the intershells
            CALL storeShells(model, intShell, totShell, siz)
            
            ! Broadcast abundances so threads are
            ! synchronized in value
            DO ii = 1, totShell - 1
                CALL MPI_BCAST(intShell(ii)%dens, siz, MPI_DOUBLE_PRECISION, &
                               0, MPI_COMM_WORLD, ierror)
            END DO
            
            ! Add the c13 extra shells
            CALL addExtraShells(intShell, extraShells, c13indx, p1indx, &
                                n14indx, totShell, extraNum, lowMass, upMass, &
                                lastConv, yscale, siz)
            
            ! Let know if C13 extension performed
            IF ((isMaster).AND.(ASSOCIATED(extraShells))) THEN
                PRINT*, "# C13 shell extension performed"
            END IF
            
            ! Create liteShell
            CALL createLiteShells(intShell, liteShell, extraShells, totShell, &
                                  nLiteShell, siz)
        END IF
        
        ! Now create the convective shell and mix
        CALL storeConvection(liteShell, nLiteShell, n1indx, siz)
        
        ! Integration dt in years
        intDt = (1.D1**t2 - 1.D1**intT1)
        
        ! Flag to know if conditions for integrating are met
        IF ((isPulse.AND.(intDt.GE.minPulsIntDt)).OR.(intDt.GE.minIntDt)) THEN
            performDt = .TRUE.
        ELSE
            performDt = .FALSE.
        END IF
        
        IF (performDt) THEN
            ! Allocate reaction array
            ALLOCATE(shellReacts(nLiteShell))
            
            ! Calculate the initial and final indices for integrating
            firstIntegShell = 0
            lastIntegShell = 0
            
            DO ii = 1, nLiteShell
                ! Flag for integrating
                integThis = .FALSE.
                
                ! Store values
                temp = liteShell(ii)%temp
                dens = liteShell(ii)%dens
                
                ! Check abundances to rule out core
                DO jj = 1, SIZE(eIndices)
                    IF (eIndices(jj)%posit.NE.0) THEN
                        jjindx = eIndices(jj)%posit
                        
                        IF (dens(jjindx).GT.yscale) THEN
                            integThis = .TRUE.
                            EXIT
                        END IF
                    END IF
                END DO
                
                IF (.NOT.integThis) CYCLE
                
                ! Exit if temperature is lower than 1e4 K or if we are in
                ! the last 10% of the star mass
                IF ((temp.LT.1.D-5).OR.(liteShell(ii)%mass0/mass.GT.0.9)) EXIT
                
                ! Update shell labels
                intgShls = intgShls + 1
                IF (firstIntegShell.EQ.0) firstIntegShell = ii
                lastIntegShell = ii
            END DO
            
            ! Calculate overshooting lower limit
            IF ((ovParam.GT.0.D0).OR.(ovPDCZParam.GT.0.D0)) THEN
                
                IF (ovMode.EQ."advective") THEN
                    CALL advOvLimit(liteShell, firstOv, p1indx, he4indx, &
                                    ovParam, ovPDCZParam)
                ELSE IF (ovMode.EQ."diffusive") THEN
                    CALL diffOvLimit(liteShell, firstOv, p1indx, he4indx, &
                                    ovParam, ovPDCZParam)
                END IF
                
                IF (firstOv.LT.firstIntegShell) firstIntegShell = firstOv
            ELSE
                firstOv = firstIntegShell
            END IF
            
            ! Allocate everything between firstIntegShell and lastIntegShell
            DO ii = firstIntegShell, lastIntegShell
                ! Divide memory allocationg among processes
                calcProc = MOD(ii - 1, nProc)
                IF (calcProc.NE.rank) CYCLE
                
                ! Store values
                temp = liteShell(ii)%temp
                rho = liteShell(ii)%rho
                
                ! Check temperature to choose reactions
                IF (temp.GE.1.D-2) THEN
                    isLowTemp = .FALSE.
                    reacts => highTempReacts
                ELSE
                    isLowTemp = .TRUE.
                    reacts => lowTempReacts
                END IF
                
                ! Create and calculate reactions
                CALL calculateReacts(shellReacts(ii), reacts, temp, rho, &
                                     partfunct, isLowTemp)
            END DO
            
            IF (isMaster) PRINT*, "# Integration timestep years: ", intDt
            
            ! dt in seconds
            intDt = intDt*sInYear
            
            ! Integrate proper
            CALL mixedIntegration(intDt, liteShell, nLiteShell, shellReacts, &
                        ntwkMass, eIndices, n1indx, p1indx, he4indx, ovParam, &
                        ovPDCZParam, mixFreq, eps, siz, yscale, firstOv, &
                        firstIntegShell, lastIntegShell, ovMode, nProc, rank)
            
            ! Reset intT1:
            intT1 = t2
            
            ! Update nSteps and print
            nSteps = nSteps + 1
            IF (isMaster) PRINT*, "# Writing in", writeFreq - nSteps, "steps"
            
            ! Deallocate reaction array
            DO ii = firstIntegShell, lastIntegShell
                ! Will only continue if correct thread
                calcProc = MOD(ii - 1, nProc)
                IF (calcProc.NE.rank) CYCLE
                
                ! Deallocate proper
                CALL clearCrossNodes(shellReacts(ii)%next)
            END DO
            DEALLOCATE(shellReacts)
        END IF
        
        ! Final mix and convection cleaning
        CALL cleanConvection(liteShell, nLiteShell, n1indx, siz)
        ! ====================================
        
        ! Write data
        IF (nSteps.GE.writeFreq) THEN
            nSteps = 0
            IF (isMaster) THEN
                ! Undo shells transformations
                CALL undoLiteShell(intShell, liteShell, siz)
                CALL undoTransformations(intShell, model, totShell)
                
                PRINT*, "# Writing model", modNum
                CALL writeData(model, totShell, modNum, mass, t2)
            END IF
        END IF
        
        ! Update t1 for the modDt calculation
        t1 = t2
        
        ! If did not integrate, keep physics
        IF (performDt) THEN
            ! Undo shells transformations
            CALL undoLiteShell(intShell, liteShell, siz)
            CALL undoTransformations(intShell, model, totShell)
            
            ! Store extraShells if they existed
            IF (ASSOCIATED(extraShells)) THEN
                ! Look for first shell
                DO jj = 1, totShell - 1
                    IF (intShell(jj)%mass0.EQ.extraShells(1)%mass0) EXIT
                END DO
                
                ! Add them
                DO ii = 1, SIZE(extraShells)
                    extraShells(ii) = intShell(jj)
                    jj = jj + 1
                END DO
            END IF
            
            ! Deallocate liteShell
            DO ii = 1, nLiteShell
                DEALLOCATE(liteShell(ii)%dens)
            END DO
            DEALLOCATE(liteShell)
            
            ! Deallocate intShell
            DO ii = 1, totShell - 1
                DEALLOCATE(intShell(ii)%dens)
            END DO
            DEALLOCATE(intShell)
            
            ! Input chemistry in new model
            CALL chemInPhys(model, nextModel, totShell, siz)
            
            ! Deallocate model
            DO ii = 1, totShell
                DEALLOCATE(model(ii)%dens)
            END DO
            DEALLOCATE(model)
            
            model => nextModel
            modNum = newModNum
            
            ! Total number of shells
            totShell = SIZE(model)
        ELSE
            ! Deallocate unused new model
            DO ii = 1, SIZE(nextModel)
                DEALLOCATE(nextModel(ii)%dens)
            END DO
            DEALLOCATE(nextModel)
        END IF
        
        CALL CPU_TIME(time2)
        
        tTot = tTot + (time2 - time1)
        IF (isMaster) THEN
            PRINT*, "# Elapsed time seconds:", time2 - time1
            
            IF (intgShls.NE.0) THEN
                PRINT*, "# Number of shells integrated:", intgShls
                PRINT*, "# Seconds per shell:", (time2 - time1)/intgShls
            END IF
            
            PRINT*, "# Total elapsed time hours:", tTot/3600.D0
            PRINT*
        END IF
    END DO
    
    ! Close physical file
    CLOSE(UNIT = uni)
    
    ! Delete continuation file
    OPEN(UNIT = uni, FILE = contChem)
    CLOSE(UNIT = uni, STATUS = "DELETE")
    
    ! Free structure arrays
    DO ii = 1, SIZE(model)
        IF (ALLOCATED(model(ii)%dens)) DEALLOCATE(model(ii)%dens)
    END DO
    
    ! Free allocatables
    DEALLOCATE(model, partfunct, dens, ntwkMass)
    
    ! Free reacts linked list
    CALL clearNodes(highTempReacts)
    CALL clearNodes(lowTempReacts)

    CALL MPI_FINALIZE(ierror)
END PROGRAM snuppat

MODULE writer_reader
    USE MPI
    USE readvars_mod
    USE structures_mod
    USE shell_operations
    IMPLICIT NONE
    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine reads the initial ATON data and puts them in model.      !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -model, SHELL type array to store the physics and chemistry.             !!!
!!! -age, first model age.                                                   !!!
!!! -localNtwk, relation between indices and names for our network.          !!!
!!! -modelNum, initial model number.                                         !!!
!!! -Metallicity Z in mass fraction, initial model number.                   !!!
!!! -n1indx, neutrons index.                                                 !!!
!!! -isPulse, holds true if in a pulse.                                      !!!
!!! -siz, the number of species (size of dens).                              !!!
!!! -rank, this thread index.                                                !!!
!!!                                                                          !!!
!!! At the output, model, age, localNtwk, and modelNum are given values.     !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialModel(model, age, localNtwk, modelNum, metallicity, n1indx, &
                        isPulse, siz, cont, rank)
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::ierror
    
    ! Input
    TYPE (SHELL), POINTER::model(:)
    DOUBLE PRECISION::age, metallicity
    INTEGER::localNtwk(:), modelNum, n1indx, siz, rank
    LOGICAL::isPulse, cont
    
    ! Local
    TYPE (SHELL), POINTER::abundncs(:)
    DOUBLE PRECISION, ALLOCATABLE::tempAbund(:), solarAbundances(:)
    INTEGER, ALLOCATABLE::relation(:), atonToNtwk(:), solarWeights(:)
    INTEGER, ALLOCATABLE::solarToNtwk(:)
    CHARACTER(3), ALLOCATABLE::names(:), names2(:), names3(:)
    CHARACTER::answer
    DOUBLE PRECISION::mass, factor, zz
    INTEGER::chemshells, atonsiz, solarsiz, ii, jj, inModNum
    LOGICAL::exst, isPhysMod, isMaster
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Identify master
    isMaster = (rank.EQ.0)
    
    ! First read chemistry, then physics and finally merge them
    
    ! Check if continuation file
    cont = .FALSE.
    IF (isMaster) THEN
        INQUIRE(FILE = contChem, EXIST = exst)
        IF (exst) THEN
            WRITE(*, '(A)', ADVANCE = "NO") "## There is a continuation file. "
            WRITE(*, '(A)', ADVANCE = "NO") "Do you want to start from its"
            WRITE(*, '(A)', ADVANCE = "NO") " model? (y/n): "
            READ(*, *) answer
            
            IF (answer.NE."n") THEN
                cont = .TRUE.
            ELSE
                ! Erase contChem
                OPEN(UNIT = uni2, FILE = contChem)
                CLOSE(UNIT = uni2, STATUS = "DELETE")
            END IF
        END IF
    END IF
    CALL MPI_BCAST(cont, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierror)
    
    ! Start from the beginning
    IF (.NOT.cont) THEN
        ! Erase previous output
        IF (isMaster) THEN
            OPEN(UNIT = uni2, FILE = output)
            CLOSE(UNIT = uni2, STATUS = "DELETE")
        END IF
        
        ! Open chemistry file, read first line and allocate everything
        OPEN(UNIT = uni2, FILE = chemistry)
        
        READ(uni2, *) inModNum, mass, chemshells
        
        ! Open aton species file to get the size
        atonsiz = 0
        OPEN(UNIT = uni3, FILE = atonSpecies)
        DO
            READ(uni3, *, IOSTAT = error)
            IF (error.NE.0) EXIT
            
            atonsiz = atonsiz + 1
        END DO
        CLOSE(UNIT = uni3)
        
        ALLOCATE(abundncs(chemshells), relation(atonsiz), tempAbund(atonsiz))
        ALLOCATE(atonToNtwk(siz), names(atonsiz), names2(siz))
        ALLOCATE(solarToNtwk(siz))
        
        ! Array relation holds the atom mass in g/mole in the position
        ! for each isotope in ATON
        OPEN(UNIT = uni3, FILE = atonSpecies)
        DO ii = 1, atonsiz
            READ(uni3, *) relation(ii), names(ii)
        END DO
        CLOSE(UNIT = uni3)
        
        ! Read the local network
        OPEN(UNIT = uni3, FILE = species)
        DO ii = 1, siz
            READ(uni3, *) localNtwk(ii), names2(ii)
        END DO
        CLOSE(UNIT = uni3)
        
        ! Read the solar abundances
        OPEN(UNIT = uni3, FILE = solar)
        
        READ(uni3, *) solarsiz
        ALLOCATE(solarAbundances(solarsiz), solarWeights(solarsiz))
        ALLOCATE(names3(solarsiz))
        
        DO ii = 1, solarsiz
            READ(uni3, *) solarWeights(ii), names3(ii), solarAbundances(ii)
        END DO
        
        CLOSE(UNIT = uni3)
        
        ! Array atonToNtwk holds the indices where Aton elements are in
        ! relation with our data
        atonToNtwk = 0
        DO ii = 1, atonsiz
            ! At each step we look for the localNtwk equivalent element and in
            ! that position in atonToNtwk we input the value ii, which is its
            ! position in relation
            DO jj = 1, siz
                IF (relation(ii).EQ.localNtwk(jj)) THEN
                    IF (names(ii).EQ.names2(jj)) THEN
                        atonToNtwk(jj) = ii
                        EXIT
                    END IF
                END IF
            END DO
        END DO
        
        ! Array solarToNtwk holds the indices where solar abundances are in
        ! relation with our data
        solarToNtwk = 0
        DO ii = 1, solarsiz
            ! At each step we look for the localNtwk equivalent element and in
            ! that position in solarToNtwk we input the value ii, which is its
            ! position in relation
            DO jj = 1, siz
                IF (solarWeights(ii).EQ.localNtwk(jj)) THEN
                    IF (names3(ii).EQ.names2(jj)) THEN
                        solarToNtwk(jj) = ii
                        EXIT
                    END IF
                END IF
            END DO
        END DO
        
        ! Fill the array
        DO ii = 1, chemshells
            ! Read line
            READ(uni2, *) abundncs(ii)%mass, tempAbund
            abundncs(ii)%mass = abundncs(ii)%mass*mass
            
            ! Convert abundances in mass fraction to abundances
            ! in number fraction (mass fraction/mass (mole/g)).
            tempAbund = tempAbund/relation
            
            ! Input ATON values.
            ALLOCATE(abundncs(ii)%dens(siz)); abundncs(ii)%dens = 0.D0
            DO jj = 1, siz
                IF (atonToNtwk(jj).NE.0) THEN
                    abundncs(ii)%dens(jj) = tempAbund(atonToNtwk(jj))
                END IF
            END DO
        END DO
        
        ! Close chemistry file
        CLOSE(UNIT = uni2)
        
        ! Now input solar values
        ! Scale abundances
        zz = 0
        DO ii = 3, solarsiz
            zz = zz + solarAbundances(ii)*solarWeights(ii)
        END DO
        factor = metallicity/zz
        
        ! CAUTION We scale all the metals the same, supposing that the key ones
        ! from the SDU come from the evolutionary model.
        solarAbundances(3:) = solarAbundances(3:)*factor
        
        DO ii =  1, chemshells
            
            ! Introduce values and keep neutrons to 0.
            DO jj = 1, siz
                IF ((atonToNtwk(jj).EQ.0).AND.(solarToNtwk(jj).NE.0)) THEN
                    abundncs(ii)%dens(jj) = solarAbundances(solarToNtwk(jj))
                END IF
            END DO
            abundncs(ii)%dens(n1indx) = 0.D0
        END DO
        
        ! Deallocate the ones used for this part
        DEALLOCATE(relation, tempAbund, atonToNtwk, names, names2)
        DEALLOCATE(solarToNtwk, solarAbundances, solarWeights, names3)
    ELSE
    ! Start from contChem.dat
        OPEN(UNIT = uni2, FILE = contChem)
        
        READ(uni2, *) inModNum, mass, chemshells
        
        ALLOCATE(abundncs(chemshells))
        
        ! Read the local network
        OPEN(UNIT = uni3, FILE = species)
        DO ii = 1, siz
            READ(uni3, *) localNtwk(ii)
        END DO
        CLOSE(UNIT = uni3)
        
        ! Fill the array
        DO ii = 1, chemshells
            ! Allocate array
            ALLOCATE(abundncs(ii)%dens(siz))
            
            ! Read line
            READ(uni2, *) abundncs(ii)%mass, abundncs(ii)%dens
            abundncs(ii)%mass = abundncs(ii)%mass*mass
        END DO
        
        ! Close chemistry file
        CLOSE(UNIT = uni2)
    END IF
    
    ! Chemistry done, now to physics
    
    ! In physics file, read first line and allocate model
    isPhysMod = readPhysMod(model, mass, age, modelNum, isPulse, inModNum, siz)
    IF (.NOT.isPhysMod) THEN
        PRINT*, "No physical information"
        STOP
    END IF
    
    ! Now introduce the abundances in model interpolating in mass
    CALL chemInPhys(abundncs, model, chemshells, siz)
    
    ! Free memory
    DO ii = 1, chemshells
        DEALLOCATE(abundncs(ii)%dens)
    END DO
    
    ! Deallocate the common ones
    DEALLOCATE(abundncs)
    
END SUBROUTINE initialModel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function reads the next data from physics. Returns true if it has   !!!
!!! been successfull and false if it hasn't.                                 !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -model, SHELL type array to store the physics and chemistry.             !!!
!!! -mass, model mass.                                                       !!!
!!! -age, model age.                                                         !!!
!!! -modelNum, model number.                                                 !!!
!!! -isPulse, holds true if in a pulse.                                      !!!
!!! -inModNum, model number to extract. If 0, extract first.                 !!!
!!! -siz, the number of species (size of dens).                              !!!
!!!                                                                          !!!
!!! At the output, everything but inModNum and siz is given the correct      !!!
!!! values.                                                                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION readPhysMod(model, mass, age, modelNum, isPulse, inModNum, siz)
    IMPLICIT NONE
    
    ! Input
    TYPE (SHELL), POINTER::model(:)
    DOUBLE PRECISION::mass, age
    INTEGER::modelNum, inModNum, siz
    LOGICAL::isPulse
    
    ! Function
    LOGICAL::readPhysMod
    
    ! Local
    TYPE (SHELL), POINTER::tempModel(:)
    DOUBLE PRECISION::rmass, rtemp, rrho, rradiat, rradius, rhp, rpres, rvel
    DOUBLE PRECISION::prevMass, repsNuc, r3AepsNuc, integNuc, integ3A
    INTEGER, ALLOCATABLE::addedIndices(:)
    INTEGER::ii, nShells, error, realNShells, iindx
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    isPulse = .FALSE.
    integNuc = 0.D0; integ3A = 0.D0
    DO
        READ(uni, *, IOSTAT = error) modelNum, mass, age, nShells
        
        IF (error.NE.0) THEN
            readPhysMod = .FALSE.
            RETURN
        END IF
        
        IF (inModNum.NE.0) THEN
            IF (inModNum.EQ.modelNum) inModNum = 0
        END IF
        
        IF (inModNum.EQ.0) ALLOCATE(tempModel(nShells), addedIndices(nShells))
        
        ! Initialize
        prevMass = -1.D2
        realNShells = 0
        
        ! Fill the array
        DO ii = 1, nShells
            IF (inModNum.EQ.0) THEN
                ALLOCATE(tempModel(ii)%dens(siz))
                READ(uni, *) rmass, rtemp, rrho, rradiat, rradius, rhp, rpres, &
                             rvel, repsNuc, r3AepsNuc
            ELSE
                READ(uni, *)
                CYCLE
            END IF
            
            ! Check if this is really a different shell
            IF (ABS(rmass - prevMass).LT.1.D-50) CYCLE
            
            ! Advance counters
            realNShells = realNShells + 1
            
            ! Count the really added indices so we can fix later the input by
            ! ignoring the zero-sized shells.
            addedIndices(realNShells) = ii
            
            ! Insert values
            tempModel(ii)%mass = rmass*mass
            tempModel(ii)%temp = 1.D1**(rtemp - 9) ! in 1e9 K
            tempModel(ii)%rho = 1.D1**rrho
            tempModel(ii)%radiat = rradiat
            tempModel(ii)%radius = rradius
            tempModel(ii)%hp = rhp
            tempModel(ii)%pressure = 1.D1**(rpres)
            tempModel(ii)%velocity = rvel
            
            ! Integrate the luminosities
            IF (prevMass.GT.0.D0) THEN
                integNuc = integNuc + repsNuc*(rmass - prevMass)
                integ3A = integ3A + r3AepsNuc*(rmass - prevMass)
            END IF
            
            prevMass = rmass
        END DO
        
        IF (inModNum.EQ.0) THEN
            ! Check if in a pulse
            IF ((integ3A/integNuc).GT.0.1) isPulse = .TRUE.
            
            ! Now fix the model
            ALLOCATE(model(realNShells))
            
            DO ii = 1, realNShells
                ! Get tempModel index and allocate dens
                iindx = addedIndices(ii)
                ALLOCATE(model(ii)%dens(siz))
                
                ! Copy values
                model(ii)%mass = tempModel(iindx)%mass
                model(ii)%temp = tempModel(iindx)%temp
                model(ii)%rho = tempModel(iindx)%rho
                model(ii)%radiat = tempModel(iindx)%radiat
                model(ii)%radius = tempModel(iindx)%radius
                model(ii)%hp = tempModel(iindx)%hp
                model(ii)%pressure = tempModel(iindx)%pressure
                model(ii)%velocity = tempModel(iindx)%velocity
            END DO
            
            ! Now free memory
            DO ii = 1, nShells
                DEALLOCATE(tempModel(ii)%dens)
            END DO
            DEALLOCATE(tempModel)
            
            readPhysMod = .TRUE.
            RETURN
        END IF
    END DO
    
END FUNCTION readPhysMod
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! The subroutine writes the data to the output and continuation files.     !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -model, SHELL type array storing the physics and chemistry.              !!!
!!! -totShell, size of model.                                                !!!
!!! -modelNum, model number.                                                 !!!
!!! -mass, model mass.                                                       !!!
!!! -age, model age.                                                         !!!
!!!                                                                          !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE writeData(model, totShell, modelNum, mass, age)
    IMPLICIT NONE
    
    ! Input
    TYPE (SHELL), POINTER::model(:)
    DOUBLE PRECISION::mass, age
    INTEGER::totShell, modelNum
    
    ! Local
    INTEGER::ii

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Write to output
    
    ! Open file
    OPEN(UNIT = uni2, FILE = output, POSITION = "APPEND", &
        FORM = "UNFORMATTED", ACCESS = "STREAM")
    
    ! First, write model number and mass
    WRITE(uni2) modelNum, mass, age, totShell, 4+SIZE(model(1)%dens)
    
    ! Now write everything else
    DO ii = 1, totShell
        WRITE(uni2) model(ii)%mass/mass, model(ii)%temp, model(ii)%rho, &
                       model(ii)%radiat, model(ii)%dens
    END DO
    
    ! Close file
    CLOSE(uni2)
    
    ! Now write to contChem
    
    ! Open file
    OPEN(UNIT = uni2, FILE = contChem)

    ! Write model, mass and number of shells
    WRITE(uni2, *) modelNum, mass, totShell
    
    ! Now write mass coordinate and abundances
    DO ii = 1, totShell
        WRITE(uni2, *) model(ii)%mass/mass, model(ii)%dens
    END DO
    
    ! Close file
    CLOSE(uni2)
    
END SUBROUTINE writeData

END MODULE writer_reader

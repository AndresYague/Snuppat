MODULE loader
    USE readvars_mod
    USE structures_mod
    USE math_routines
    IMPLICIT NONE
    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the cross sections linked list.                  !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -highTempReacts, a pointer to type REACT.                                !!!
!!! -lowTempReacts, a pointer to type REACT.                                 !!!
!!!                                                                          !!!
!!! At the output, both highTempReacts and lowTempReacts are the heads of    !!!
!!! two linked lists, the first one containing all the reactions in the      !!!
!!! network, and the second one containing just decay reactions, which can   !!!
!!! be followed in low temperature shells.                                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE loadNetwork(highTempReacts, lowTempReacts)
    IMPLICIT NONE
    
    ! Input variables
    TYPE (REACT), POINTER::highTempReacts, lowTempReacts
    
    ! Local variables
    TYPE (REACT), POINTER::curr, last
    DOUBLE PRECISION, DIMENSION(7)::avector
    DOUBLE PRECISION, ALLOCATABLE::locTempTable(:)
    CHARACTER(10)::source
    INTEGER, DIMENSION(3)::neg
    INTEGER, DIMENSION(4)::pos
    INTEGER::kk, jj, jumpSiz, locTabSiz, kkindx

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! This reads the files with the reactions and makes the linked list.
    DO kk = 1, SIZE(filenames) + SIZE(lowTempFiles)
        IF ((kk - SIZE(filenames)).LT.1) THEN
            kkindx = kk
            OPEN(UNIT = uni, FILE = "data/"//filenames(kk)//".lst")
        ELSE
            kkindx = kk - SIZE(filenames)
            OPEN(UNIT = uni, FILE = "data/"//lowTempFiles(kkindx)//".lst")
        END IF
        
        ! For the first list we must skip the first jumpSiz lines.
        IF (kk.EQ.1) THEN
            READ(uni, *) jumpSiz
            
            DO jj = 1, jumpSiz
                READ(uni, *)
            END DO
        END IF
        
        ! Nullify last if we are at the head
        IF (kkindx.EQ.1) NULLIFY(last)
        
        ! In this loop the coefficients are read.
        DO
            ! This line reads the nucleon tag, the reaction source and coefs.
            neg = 0; pos = 0
            IF (kk.EQ.kkindx) THEN
                READ(uni, *, IOSTAT = error) neg(1:targ(kk)), pos(1:prod(kk)), &
                                             source, avector
            ELSE
                READ(uni, *, IOSTAT = error) neg(1:targ(kkindx)), &
                                             pos(1:prod(kkindx)), source, &
                                             locTabSiz
            END IF
            IF (error.NE.0) EXIT
            
            ! Check that tabSiz has not changed
            IF ((kk.NE.kkindx).AND.(tabSiz.NE.locTabSiz)) THEN
                PRINT*, "Not the same table size for starlib reactions!"
                PRINT*, tabSiz, locTabSiz
                STOP
            END IF
            
            ! Put reaction in the list:
            ALLOCATE(curr)
            IF (kk.EQ.kkindx) THEN
                !Allocate space for avector
                ALLOCATE(curr%avector(7))
                curr%avector = avector
            ELSE
                ! Allocate space for the tables
                ALLOCATE(locTempTable(locTabSiz))
                ALLOCATE(curr%crossTable(tabSiz))
                
                ! Read the tables
                READ(uni, *) locTempTable
                READ(uni, *) curr%crossTable
                
                ! Make all zeros in crossTable 1e-100
                DO jj = 1, tabSiz
                    IF (curr%crossTable(jj).LT.1.D-100) THEN
                        curr%crossTable(jj) = 1.D-100
                    END IF
                END DO
                
                ! Define tempTable if not done already
                IF (MAXVAL(tempTable).LE.0.D0) tempTable = locTempTable
                
                ! Check that the temperature table is the same
                IF (MAXVAL(ABS(tempTable - locTempTable)).GT.0.D0) THEN
                    PRINT*, "Not the same temperatures for starlib reactions!"
                    
                    DO jj = 1, tabSiz
                        PRINT*, tempTable(jj), locTempTable(jj)
                    END DO
                    
                    STOP
                END IF
                
                DEALLOCATE(locTempTable)
            END IF
            
            ! Input all common values
            curr%source = source
            curr%targnum = targ(kkindx)
            curr%prodnum = prod(kkindx)
            curr%totnum = targ(kkindx) + prod(kkindx)
            curr%targindx = neg
            curr%prodindx = pos
            
            IF (curr%source(1:2).EQ."ec") THEN
                curr%isEc = .TRUE.
            ELSE
                curr%isEc = .FALSE.
            END IF
            
            ! Introduce new reaction at the end
            IF (ASSOCIATED(last)) THEN
                last%next => curr
            ELSE IF (kk.EQ.1) THEN
                highTempReacts => curr
            ELSE IF (kkindx.EQ.1) THEN
                lowTempReacts => curr
            END IF
            last => curr
            
            ! Nullify next to be in the safe side
            NULLIFY(curr%next)
            
        END DO
        
        CLOSE(UNIT = uni)
    END DO
    
END SUBROUTINE loadNetwork

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine loads the partition functions in "partfunct".            !!!
!!!                                                                          !!!
!!! The input value is:                                                      !!!
!!! -partfunct, a two-dimensional array.                                     !!!
!!!                                                                          !!!
!!! At the output, the array partfunct contains all the partition function   !!!
!!! tabulated values for a given nucleon. The first index is the label of    !!!
!!! said nucleon.                                                            !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE loadPartitions(partfunct)
    IMPLICIT NONE
    
    ! Input variables
    DOUBLE PRECISION::partfunct(:, :)
    
    ! Local variables
    INTEGER::ii
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Read file and store
    OPEN(UNIT = uni, FILE = "data/partition.lst")
    DO
        READ(uni, *, IOSTAT = error) ii
        READ(uni, *, IOSTAT = error) partfunct(ii, :)
        IF (error.NE.0) EXIT
    END DO
    CLOSE(UNIT = uni)
    
END SUBROUTINE loadPartitions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the cross section for every reaction for the  !!!
!!! given temperature and density.                                           !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -reacts, the linked list of CROSSARR types with the cross sections.      !!!
!!! -fullReacts, the linked list of REACT types with the cross sections.     !!!
!!! -temp, the temperature in T9 units.                                      !!!
!!! -rho, the density in g/cm^3.                                             !!!
!!! -partfun, the array with the partition function tables.                  !!!
!!! -isLowTemp, boolean telling us if we are in low temperature reactions.   !!!
!!!                                                                          !!!
!!! At the output, all the cross sections in the linked list have been       !!!
!!! calculated for the given temperature and density.                        !!!
!!!                                                                          !!!
!!! An additional function is included to calculate the cross sections.      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculateReacts(reacts, fullReacts, temp, rho, partfun, isLowTemp)
    IMPLICIT NONE
    
    ! Input variables
    TYPE (CROSSARR), TARGET::reacts
    TYPE (REACT), TARGET::fullReacts
    DOUBLE PRECISION::temp, rho, partfun(:, :)
    LOGICAL::isLowTemp
    
    ! Local variables
    TYPE (CROSSARR), POINTER::cCross
    TYPE (REACT), POINTER::cReact
    DOUBLE PRECISION::crsect, partval
    INTEGER::ii, jj1, jj2, jj3, targIndx(3), prodIndx(4)
    LOGICAL::isRepeated, firstAdded, isSameReaction

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize
    targIndx = 0
    prodIndx = 0
    
    ! Point to lists beginning
    cCross => reacts
    cReact => fullReacts
    firstAdded = .FALSE.
    
    ! In this loop the reactions are calculated
    DO
        ! Calculate cross section.
        IF (isLowTemp) THEN
            crsect = interpolateOneValue(log(temp), log(tempTable), &
                                         log(cReact%crossTable))
            crsect = exp(crsect)
        ELSE
            crsect = highTempCross(cReact%avector, temp)
        END IF
        crsect = crsect*(rho**(cReact%targnum - 1))
        
        ! If reaction is "ec", we have to multiply again by rho.
        IF (cReact%isEc) crsect = crsect*rho
        
        ! If there are more than one of each target, then we have to divide
        ! by the factorial of the number of repetitions, so by 2 or by 6.
        IF (cReact%targnum.GE.2) THEN
            jj1 = cReact%targindx(1)
            jj2 = cReact%targindx(2)
            jj3 = cReact%targindx(3)
            
            IF ((jj1.EQ.jj2).AND.(jj2.EQ.jj3)) THEN
                crsect = crsect/6.D0
            ELSE IF ((jj1.EQ.jj2).OR.(jj1.EQ.jj3).OR.(jj2.EQ.jj3)) THEN
                crsect = crsect/2.D0
            END IF
        END IF
        
        ! This block reads the partition function and calculates it.
        IF ((cReact%source(5:5).EQ.'v').OR.(cReact%source(6:6).EQ.'v')) THEN
            
            ! If the temperature is lower than 10^8 K, the partition
            ! functions will always be 1, if the temperature is higher,
            ! a logarithmical interpolation is made.
            IF ((temp.GE.1D-1).AND.(.NOT.isLowTemp)) THEN
                DO ii = 1, cReact%targnum + cReact%prodnum
                    IF (ii.LE.cReact%targnum) THEN
                        jj1 = cReact%targindx(ii)
                        
                        CALL partitionValue(temp, partfun(jj1, :), partval)
                        crsect = crsect/partval
                    ELSE
                        jj1 = cReact%prodindx(ii - cReact%targnum)
                        
                        CALL partitionValue(temp, partfun(jj1, :), partval)
                        crsect = crsect*partval
                    END IF
                END DO
            END IF
        END IF
        
        ! Check that crsect is bigger than a value
        ! If the reaction is a repeat, simply add it to the last one
        IF (crsect.GT.1.D-40) THEN
            isRepeated = .TRUE.
            isSameReaction = .TRUE.
            IF (targIndx(1).NE.cReact%targIndx(1)) THEN
                isRepeated = .FALSE.
                isSameReaction = .FALSE.
            ELSE IF (targIndx(2).NE.cReact%targIndx(2)) THEN
                isRepeated = .FALSE.
                isSameReaction = .FALSE.
            ELSE IF (targIndx(3).NE.cReact%targIndx(3)) THEN
                isRepeated = .FALSE.
                isSameReaction = .FALSE.
            ELSE IF (prodIndx(1).NE.cReact%prodIndx(1)) THEN
                isSameReaction = .FALSE.
            ELSE IF (prodIndx(2).NE.cReact%prodIndx(2)) THEN
                isSameReaction = .FALSE.
            ELSE IF (prodIndx(3).NE.cReact%prodIndx(3)) THEN
                isSameReaction = .FALSE.
            ELSE IF (prodIndx(4).NE.cReact%prodIndx(4)) THEN
                isSameReaction = .FALSE.
            END IF
            
            ! If not repeated, fill this bit and allocate the next
            ! If repeated, add to the last one
            IF (.NOT.isSameReaction) THEN
                IF (firstAdded) THEN
                    ALLOCATE(cCross%next)
                    cCross => cCross%next
                    NULLIFY(cCross%next)
                ELSE
                    firstAdded = .TRUE.
                END IF
                
                cCross%crossect = crsect
                cCross%targIndx = cReact%targIndx
                cCross%prodIndx = cReact%prodIndx
                cCross%targnum = cReact%targnum
                cCross%prodnum = cReact%prodnum
                cCross%totnum = cReact%totnum
                cCross%isEc = cReact%isEc
                cCross%isRepeated = isRepeated
            ELSE
                cCross%crossect = cCross%crossect + crsect
            END IF
            
            targIndx = cReact%targIndx
            prodIndx = cReact%prodIndx
        END IF
        
        ! Go to next reaction
        IF (ASSOCIATED(cReact%next)) THEN
            cReact => cReact%next
        ELSE
            EXIT
        END IF
    END DO
    
END SUBROUTINE calculateReacts

FUNCTION highTempCross(avector, temp) RESULT(crsect)
    IMPLICIT NONE
    
    ! Input variables
    DOUBLE PRECISION::avector(:), temp
    
    ! Function
    DOUBLE PRECISION::crsect
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    crsect = avector(1) + avector(2)/temp + avector(3)*temp**(-1.d0/3.d0) + &
             avector(4)*temp**(1.d0/3.d0) + avector(5)*temp + &
             avector(6)*temp**(5.d0/3.d0) + avector(7)*log(temp)
    
    ! Check that there are no underflows
    IF (crsect.LT.-300) THEN
        crsect = 0
    ELSE
        crsect = exp(crsect)
    END IF
    
END FUNCTION highTempCross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the partition function for each element using !!!
!!! the data from reaclib and performing a logarithmic interpolation.        !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -temp, the temperature in T9 units.                                      !!!
!!! -partit, the array with the partition function tables.                   !!!
!!! -partVal, the interpolated value.                                        !!!
!!!                                                                          !!!
!!! At the output, partVal carries the interpolated value.                   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE partitionValue(temp, partit, partVal)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::temp, partit(24), partVal
    
    ! Temperature array for interpolation (in T9).
    DOUBLE PRECISION, DIMENSION(24)::tempArr = (/ 1D-1, 1.5D-1, 2D-1, 3D-1, &
                4D-1, 5D-1, 6D-1, 7D-1, 8D-1, 9D-1, 1D0, 1.5D0, 2D0, 2.5D0, &
                3D0, 3.5D0, 4D0, 4.5D0, 5D0, 6D0,7D0, 8D0, 9D0, 1D1 /)

    ! Local
    INTEGER::kk
    DOUBLE PRECISION::aa, bb
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! First we search between which two temperatures is our temp.
    DO kk = 1, 23
        IF ((temp.GE.tempArr(kk)).AND.(temp.LE.tempArr(kk+1))) EXIT
    END DO
    
    ! Now we avoid calculation if the partition function is constant between the
    ! two temperatures and make the logarithmic interpolation if it's not.
    IF (ABS(partit(kk) - partit(kk+1)).LE.1.D-5) THEN
        partVal = partit(kk)
    ELSE
        aa = (log(partit(kk+1)) - log(partit(kk)))/ &
             (log(tempArr(kk+1)) - log(tempArr(kk)))
        
        bb = -(log(partit(kk+1))*log(tempArr(kk)) - log(partit(kk))* &
              log(tempArr(kk+1)))/(log(tempArr(kk+1)) - log(tempArr(kk)))
        
        partVal = (temp**aa)*exp(bb)
    END IF

END SUBROUTINE partitionValue

END MODULE loader

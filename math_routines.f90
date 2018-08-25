MODULE math_routines
    USE structures_mod
    IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function interpolates using a cubic monotone piecewise              !!!
!!! interpolation. It takes the x array we want to interpolate to (mass),    !!!
!!! the known x values (valMass), array (valArr), and their size (tot).      !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -mass, the array with mass coordinates for which we ignore the y values. !!!
!!! -valMass, the array with mass coordinates with known "y" values.         !!!
!!! -valARR, the known "y" array.                                            !!!
!!! -tot, the size of the mass array.                                        !!!
!!!                                                                          !!!
!!! Returns the interpolated values.                                         !!!
!!!                                                                          !!!
!!! Ref: Steffen, M. 1990, A&A, 239, 443.                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolate(mass, valMass, valArr, tot)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::mass(:), valMass(:), valArr(:)
    INTEGER::tot
    
    ! Function
    DOUBLE PRECISION::interpolate(tot)
    
    ! Local
    DOUBLE PRECISION::yp1, y0, ym1, yp2, si0, sip1, sim1, deriv0, pip1, ai, bi
    DOUBLE PRECISION::derivp1, pi0, xp1, xp2, x0, xm1, hi0, hip1, him1, deltx
    INTEGER::p1, p2, ii, jj, valShells
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! The interpolated value is given by the formula:
    ! f_i(x) = ai*(x - x_i)^3 + bi*(x - x_i)^2 + ci*(x - x_i) + di
    
    ! Initialize some variables so the compiler does not complain
    xp1 = 0; x0 = 0; p1 = 1; p2 = 1
    
    ! Size for valMass
    valShells = SIZE(valMass)
    
    DO ii = 1, tot
        ! Get p1 and p2 for this mass and store x0 and xp1
        DO jj = p1, valShells - 1
            x0 = valMass(jj)
            xp1 = valMass(jj + 1)
            
            IF ((x0.LE.mass(ii)).AND.(xp1.GE.mass(ii))) THEN
                p1 = jj; p2 = jj + 1
                EXIT
            END IF
        END DO
        
        ! Fix p1 and p2 status if not correct
        IF (xp1.LT.mass(ii)) THEN
            p1 = valShells - 1
            p2 = valShells
        END IF
        
        ! Store yp1 and y0
        yp1 = valArr(p2); y0 = valArr(p1)
        
        hi0 = (xp1 - x0)
        si0 = (yp1 - y0)/hi0
        
        ! Calculate sim1 and deriv0
        IF (p1.GT.1) THEN
            ! Store x, y
            ym1 = valArr(p1 - 1)
            xm1 = valMass(p1 - 1)
            
            him1 = (x0 - xm1)
            sim1 = (y0 - ym1)/him1
            
            ! Pi0 calculation
            pi0 = (sim1*hi0 + si0*him1)/(him1 + hi0)
            
            ! Derivative
            deriv0 = SIGN(1.D0, sim1) + SIGN(1.D0, si0)
            deriv0 = deriv0*MIN(ABS(sim1), ABS(si0), 0.5D0*ABS(pi0))
        ELSE
            ! We are in the lowest extreme
            deriv0 = 0.D0
        END IF
        
        ! Calculate sip1, pip1 and derivp1
        IF (p2.LT.valShells) THEN
            yp2 = valArr(p2 + 1)
            xp2 = valMass(p2 + 1)
            
            hip1 = (xp2 - xp1)
            sip1 = (yp2 - yp1)/hip1
            
            ! Pip1 calculation
            pip1 = (si0*hip1 + sip1*hi0)/(hi0 + hip1)
            
            ! Derivative
            derivp1 = SIGN(1.D0, si0) + SIGN(1.D0, sip1)
            derivp1 = derivp1*MIN(ABS(si0), ABS(sip1), 0.5D0*ABS(pip1))
        ELSE
            ! We are in the highest extreme
            derivp1 = 0.D0
        END IF
        
        ! Now calculate coefficients (ci = deriv0; di = y0)
        ai = (deriv0 + derivp1 - 2.D0*si0)/(hi0*hi0)
        bi = (3.D0*si0 - 2*deriv0 - derivp1)/hi0
        
        deltx = mass(ii) - x0
        interpolate(ii) = ai*deltx**3 + bi*deltx**2 + deriv0*deltx + y0
    END DO
    
    RETURN
    
END FUNCTION interpolate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Same as interpolate but for a scalar.                                    !!!
!!!                                                                          !!!
!!! Ref: Steffen, M. 1990, A&A, 239, 443.                                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION interpolateOneValue(mass, valMass, valArr)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::mass, valMass(:), valArr(:)
    
    ! Function
    DOUBLE PRECISION::interpolateOneValue
    
    ! Local
    DOUBLE PRECISION::yp1, y0, ym1, yp2, si0, sip1, sim1, deriv0, pip1, ai, bi
    DOUBLE PRECISION::derivp1, pi0, xp1, xp2, x0, xm1, hi0, hip1, him1, deltx
    INTEGER::p1, p2, jj, valShells
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! The interpolated value is given by the formula:
    ! f_i(x) = ai*(x - x_i)^3 + bi*(x - x_i)^2 + ci*(x - x_i) + di
    
    ! Initialize some variables so the compiler does not complain
    xp1 = 0; x0 = 0; p1 = 1; p2 = 1
    
    ! Size for valMass
    valShells = SIZE(valMass)
    
    ! Get p1 and p2 for this mass and store x0 and xp1
    DO jj = 1, valShells - 1
        x0 = valMass(jj)
        xp1 = valMass(jj + 1)
        
        IF ((x0.LE.mass).AND.(xp1.GE.mass)) THEN
            p1 = jj; p2 = jj + 1
            EXIT
        END IF
    END DO
    
    ! Fix p1 and p2 status if not correct
    IF (xp1.LT.mass) THEN
        p1 = valShells - 1
        p2 = valShells
    END IF
    
    ! Store yp1 and y0
    yp1 = valArr(p2); y0 = valArr(p1)
    
    hi0 = (xp1 - x0)
    si0 = (yp1 - y0)/hi0
    
    ! Calculate sim1 and deriv0
    IF (p1.GT.1) THEN
        ! Store x, y
        ym1 = valArr(p1 - 1)
        xm1 = valMass(p1 - 1)
        
        him1 = (x0 - xm1)
        sim1 = (y0 - ym1)/him1
        
        ! Pi0 calculation
        pi0 = (sim1*hi0 + si0*him1)/(him1 + hi0)
        
        ! Derivative
        deriv0 = SIGN(1.D0, sim1) + SIGN(1.D0, si0)
        deriv0 = deriv0*MIN(ABS(sim1), ABS(si0), 0.5D0*ABS(pi0))
    ELSE
        ! We are in the lowest extreme
        deriv0 = 0.D0
    END IF
    
    ! Calculate sip1, pip1 and derivp1
    IF (p2.LT.valShells) THEN
        yp2 = valArr(p2 + 1)
        xp2 = valMass(p2 + 1)
        
        hip1 = (xp2 - xp1)
        sip1 = (yp2 - yp1)/hip1
        
        ! Pip1 calculation
        pip1 = (si0*hip1 + sip1*hi0)/(hi0 + hip1)
        
        ! Derivative
        derivp1 = SIGN(1.D0, si0) + SIGN(1.D0, sip1)
        derivp1 = derivp1*MIN(ABS(si0), ABS(sip1), 0.5D0*ABS(pip1))
    ELSE
        ! We are in the highest extreme
        derivp1 = 0.D0
    END IF
    
    ! Now calculate coefficients (ci = deriv0; di = y0)
    ai = (deriv0 + derivp1 - 2.D0*si0)/(hi0*hi0)
    bi = (3.D0*si0 - 2*deriv0 - derivp1)/hi0
    
    deltx = mass - x0
    interpolateOneValue = ai*deltx**3 + bi*deltx**2 + deriv0*deltx + y0
    
    RETURN
    
END FUNCTION interpolateOneValue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the jacobian.                                    !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -crosLst, a linked list of CROSSARR types with the network cross         !!!
!!!            sections.                                                     !!!
!!! -dens, the shell abundances.                                             !!!
!!! -jacobian, the output jacobian matrix.                                   !!!
!!! -eIndices, an array of type INDICES with the indices of the maximum      !!!
!!!            contributors to the free electrons in the star.               !!!
!!!                                                                          !!!
!!! On output the jacobian is filled with the derivatives of the rates.      !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createJacobian(crosLst, dens, jacobian, eIndices)
    IMPLICIT NONE

    ! Input
    TYPE (CROSSARR), TARGET::crosLst
    TYPE (INDICES)::eIndices(:)
    DOUBLE PRECISION::dens(:), jacobian(:, :)
    
    ! Local
    TYPE (CROSSARR), POINTER::cCross
    DOUBLE PRECISION::derivative, Ye
    INTEGER::ii, jj, jjindx, iiindx
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    cCross => crosLst
    DO
        Ye = 1.D0
        
        ! If the reaction is "ec" we have to calculate Ye.
        IF (cCross%isEc) CALL calculateYe(dens, eIndices, Ye)
        
        ! Now calculate the jacobian elements: jacobian(j, i) = df_j/dy_i
        ! Bear in mind we are creating it transposed
        DO ii = 1, cCross%targnum
            iiindx = cCross%targindx(ii)
            
            ! First calculate the derivative. We have to multiply
            ! by the other targets.
            derivative = cCross%crossect*Ye
            DO jj = 1, cCross%targnum
                IF (jj.NE.ii) THEN
                    jjindx = cCross%targindx(jj)
                    derivative = derivative*dens(jjindx)
                END IF
            END DO
         
            ! Now input it in the jacobian.
            DO jj = 1, cCross%targnum + cCross%prodnum
                IF (jj.LE.cCross%targnum) THEN
                    jjindx = cCross%targindx(jj)
                    jacobian(iiindx, jjindx) = jacobian(iiindx, jjindx) - &
                                               derivative
                ELSE
                    jjindx = cCross%prodindx(jj - cCross%targnum)
                    jacobian(iiindx, jjindx) = jacobian(iiindx, jjindx) + &
                                               derivative
                END IF
            END DO
        END DO
        
        IF (ASSOCIATED(cCross%next)) THEN
            cCross => cCross%next
        ELSE
            EXIT
        END IF
    END DO
    
END SUBROUTINE createJacobian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the rates array.                                 !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -crosLst, a linked list of CROSSARR types with the network cross         !!!
!!!            sections.                                                     !!!
!!! -dens, the shell abundances.                                             !!!
!!! -rates, the output rates array.                                          !!!
!!! -eIndices, an array of type INDICES with the indices of the maximum      !!!
!!!             contributors to the free electrons in the star.              !!!
!!!                                                                          !!!
!!! On output the rates array is filled with the rates values.               !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createRates(crosLst, dens, kk, dd, eIndices)
    IMPLICIT NONE

    ! Input
    TYPE (CROSSARR), TARGET::crosLst
    TYPE (INDICES)::eIndices(:)
    DOUBLE PRECISION::dens(:), kk(:), dd(:), ymin = 1.D-250
    
    ! Local
    TYPE (CROSSARR), POINTER::cCross
    DOUBLE PRECISION::rate, Ye, prod, oldProd, invDens(SIZE(dens))
    INTEGER::ii, iiindx
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Zero kk and dd
    kk = 0.D0; dd = 0.D0
    invDens = 1/(dens + ymin)
    
    ! Link first reaction
    cCross => crosLst
    
    oldProd = 1
    DO
        ! If the reaction is "ec" we have to calculate Ye.
        IF (cCross%isEc) THEN
            Ye = 1.D0
            CALL calculateYe(dens, eIndices, Ye)
            rate = cCross%crossect*Ye
        ELSE
            rate = cCross%crossect
        END IF
        
        ! Calculate the rate: rates(ii) = f(ii)
        IF (cCross%isRepeated) THEN
            prod = oldProd
        ELSE
            prod = 1
            DO ii = 1, cCross%targnum
                iiindx = cCross%targindx(ii)
                prod = prod*(dens(iiindx) + ymin)
            END DO
            oldProd = prod
        END IF
        rate = rate*prod
        
        ! Input the rate.
        DO ii = 1, cCross%targnum
            iiindx = cCross%targindx(ii)
            dd(iiindx) = dd(iiindx) + rate
        END DO
        DO ii = 1, cCross%prodnum
            iiindx = cCross%prodindx(ii)
            kk(iiindx) = kk(iiindx) + rate
        END DO
        
        IF (ASSOCIATED(cCross%next)) THEN
            cCross => cCross%next
        ELSE
            EXIT
        END IF
    END DO
    
    dd = dd*invDens
    
END SUBROUTINE createRates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates Ye.                                           !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -dens, the shell abundances.                                             !!!
!!! -eIndx, an array of type INDICES with the indices of the maximum         !!!
!!!         contributors to the free electrons in the star.                  !!!
!!! -Ye, the free electrons in mass number.                                  !!!
!!!                                                                          !!!
!!! On output, Ye has a value between 0.5 and 1 in [mole/g].                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculateYe(dens, eIndx, Ye)
    IMPLICIT NONE
    
    ! Input
    TYPE (INDICES)::eIndx(:)
    DOUBLE PRECISION::dens(:), Ye
    
    ! Local
    TYPE (INDICES)::relation(5)
    DOUBLE PRECISION::metalDens
    INTEGER::ii, jj, iiindx, zz
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Fill "relation" array, were posit is the atomic weight:
    relation(1)%elem = "p"; relation(1)%posit = 1
    relation(2)%elem = "d"; relation(2)%posit = 2
    relation(3)%elem = "t"; relation(3)%posit = 3
    relation(4)%elem = "he3"; relation(4)%posit = 3
    relation(5)%elem = "he4"; relation(5)%posit = 4
    
    ! Calculate metalicity (Xz)
    ! Here zz is A_i for p, d, he3 and he4 (zz = 1, 2, 3, 4)
    zz = 0                   ! (so the compiler doesn't complain)
    metalDens = 1.D0
    DO ii = 1, SIZE(eIndx)
        iiindx = eIndx(ii)%posit
        IF (iiindx.NE.0) THEN
            DO jj = 1, SIZE(relation)
                IF (eIndx(ii)%elem.EQ.relation(jj)%elem) THEN
                    zz = relation(jj)%posit
                    EXIT
                END IF
            END DO
            
            metalDens = metalDens - zz*dens(iiindx)
        END IF
    END DO
    
    ! Calculate Ye
    ! The electron mass number because of a metal is half the metal
    ! mass fraction (approximately)
    Ye = metalDens*5.D-1
    DO ii = 1, SIZE(eIndx)
        iiindx = eIndx(ii)%posit
        IF (iiindx.NE.0) THEN
            IF ((eIndx(ii)%elem.NE."he3").AND.(eIndx(ii)%elem.NE."he4")) THEN
                
                Ye = Ye + dens(iiindx)
            ELSE
                Ye = Ye + 2*dens(iiindx)
            END IF
        END IF
    END DO
    
END SUBROUTINE calculateYe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the system index.                                !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -matrix, the matrix for the system index.                                !!!
!!! -sysInd, the integer array with indices.                                 !!!
!!! -he4indx, index for he4.                                                 !!!
!!! -siz, the matrix dimension.                                              !!!
!!!                                                                          !!!
!!! On output, sysInd has been filled with the nonzero indices.              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createSysInd(matrix, sysInd, he4indx, siz)
    IMPLICIT NONE

    ! Input
    DOUBLE PRECISION::matrix(:, :)
    INTEGER::sysInd(:, :), he4indx, siz
    
    ! Local
    INTEGER::jj, kk, k1, k2, k3, k4, currStreak, maxStreak
    LOGICAL::inStreak
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Store nonzero indices remember is already transposed
    DO kk = 1, siz
        ! Initialize
        sysInd(kk, :) = kk
        
        ! Introduce from helium onwards
        DO jj = he4indx + 1, kk - 1
            IF (ABS(matrix(jj, kk)).GT.0.D0) THEN
                sysInd(kk, 1) = jj
                EXIT
            END IF
        END DO
        
        ! And from the end, backwards
        DO jj = siz, kk + 1, -1
            IF (ABS(matrix(jj, kk)).GT.0.D0) THEN
                sysInd(kk, 4) = jj
                EXIT
            END IF
        END DO
    END DO
    
    ! Store biggest zero streak
    DO kk = 1, siz
        k2 = 0; k3 = 0
        currStreak = 0
        sysInd(kk, 2:3) = 0
        k1 = sysInd(kk, 1); k4 = sysInd(kk, 4)
        
        inStreak = .FALSE.
        maxStreak = 0
        DO jj = k1, k4
            ! Check for maximum streak and store indices if found
            IF ((ABS(matrix(jj, kk)).GT.0.D0).OR.(jj.EQ.kk)) THEN
                IF (inStreak) THEN
                    IF (currStreak.GT.maxStreak) THEN
                        maxStreak = currStreak
                        sysInd(kk, 2) = k2
                        sysInd(kk, 3) = k3
                    END IF
                    
                    inStreak = .FALSE.
                END IF
            ELSE
                ! The indices are those of non-zero numbers
                IF (.NOT.inStreak) THEN
                    k2 = jj - 1
                    currStreak = 0
                END IF
                k3 = jj + 1
                
                currStreak = currStreak + 1
                inStreak = .TRUE.
            END IF
        END DO
        
        ! If k2 and k3 are 0, fix to avoid ifs
        k2 = sysInd(kk, 2)
        IF (k2.EQ.0) THEN
            sysInd(kk, 2:3) = k4
            sysInd(kk, 4) = k1 - 1
        END IF
    END DO
    
END SUBROUTINE createSysInd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine performs a sparse gaussian elimination.                  !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -aa, the matrix in which perform the elimination.                        !!!
!!! -tMatx, the matrix holding all the transformations done during the       !!!
!!!         elimination. That way we can transform the independent vector.   !!!
!!! -sysInd, array holding the system consecutive nonzero indices.           !!!
!!! -tInd, array holding the tMatx consecutive nonzero indices.              !!!
!!! -he4indx, index at which he4 is.                                         !!!
!!! -siz, aa size.                                                           !!!
!!!                                                                          !!!
!!! On output, aa, tMatx, sysInd and tInd have been modified.                !!!
!!!                                                                          !!!
!!! This subroutine is designed for a TRANSPOSED matrix. Transposed outside. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sparseLU(aa, tMatx, sysInd, tInd, he4indx, siz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::aa(:, :), tMatx(:, :)
    INTEGER::sysInd(:, :), tInd(:), he4indx, siz
    
    ! Local
    DOUBLE PRECISION::div, multp
    INTEGER::ii, jj, i1, i2, i3, i4
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Add index to tInd
    DO ii = 1, siz
        tInd(ii) = ii
    END DO
    
    ! Perform factorization
    DO ii = siz, 1, -1
        ! Get indices for this row
        i1 = sysInd(ii, 1); i2 = sysInd(ii, 2)
        i3 = sysInd(ii, 3); i4 = sysInd(ii, 4)
        
        ! Only operate up to the diagonal
        IF (i2.GT.ii) i2 = ii
        IF (i4.GT.ii) i4 = ii
        
        ! First, divide this row by its diagonal value
        DO WHILE(ABS(1.D0 - aa(ii, ii)).GE.1.D-16)
            div = 1.D0/aa(ii, ii)
            
            IF (ii.LE.he4indx) THEN
                aa(1:i2, ii) = aa(1:i2, ii)*div
                aa(i3:i4, ii) = aa(i3:i4, ii)*div
            ELSE
                aa(1:he4indx, ii) = aa(1:he4indx, ii)*div
                aa(i1:i2, ii) = aa(i1:i2, ii)*div
                aa(i3:i4, ii) = aa(i3:i4, ii)*div
            END IF
            
            tMatx(ii:tInd(ii), ii) = tMatx(ii:tInd(ii), ii)*div
        END DO
        
        ! Now eliminate the non diagonal values of column ii
        DO jj = ii - 1, 1, -1
            IF ((sysInd(jj, 4).GE.ii).OR.(sysInd(jj, 2).GE.ii)) THEN
                IF (ABS(aa(ii, jj)).LT.1.D-300) CYCLE
                
                ! Operation on aa
                multp = aa(ii, jj)
                IF (ii.LE.he4indx) THEN
                    aa(1:i2, jj) = aa(1:i2, jj) - aa(1:i2, ii)*multp
                    aa(i3:i4, jj) = aa(i3:i4, jj) - aa(i3:i4, ii)*multp
                ELSE
                    aa(1:he4indx, jj) = aa(1:he4indx, jj) - &
                                        aa(1:he4indx, ii)*multp
                    aa(i1:i2, jj) = aa(i1:i2, jj) - aa(i1:i2, ii)*multp
                    aa(i3:i4, jj) = aa(i3:i4, jj) - aa(i3:i4, ii)*multp
                END IF
                
                ! Operation on tMatx
                tMatx(ii:tInd(ii), jj) = tMatx(ii:tInd(ii), jj) - &
                                         tMatx(ii:tInd(ii), ii)*multp
                
                ! Change aa indices
                sysInd(jj, 1) = MIN(sysInd(ii, 1), sysInd(jj, 1))
                sysInd(jj, 2) = MAX(sysInd(ii, 2), sysInd(jj, 2))
                sysInd(jj, 3) = MIN(sysInd(ii, 3), sysInd(jj, 3))
                sysInd(jj, 4) = MAX(sysInd(ii, 4), sysInd(jj, 4))
                
                ! Fix overlaps
                IF (sysInd(jj, 2).GE.sysInd(jj, 3)) THEN
                    IF (sysInd(jj, 4).GE.sysInd(jj, 3)) THEN
                        sysInd(jj, 2) = sysInd(jj, 3) - 1
                    END IF
                END IF
                
                ! Change tMatx indices
                tInd(jj) = MAX(tInd(ii), tInd(jj))
            END IF
        END DO
    END DO
    
END SUBROUTINE sparseLU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine performs a sparse gaussian elimination tailored for the  !!!
!!! downwards overshooting matrix.                                           !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -aa, the matrix in which perform the elimination.                        !!!
!!! -tMatx, the matrix holding all the transformations done during the       !!!
!!!         elimination. That way we can transform the independent vector.   !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the furthermost index they affect.                 !!!
!!! -siz, aa size.                                                           !!!
!!!                                                                          !!!
!!! On output, aa and tMatx have been modified.                              !!!
!!!                                                                          !!!
!!! This subroutine is designed for a TRANSPOSED matrix. Transposed outside. !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sparseOvMatrixGauss(aa, tMatx, convecIndex, siz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::aa(:, :), tMatx(:, :)
    INTEGER::convecIndex(:, :)
    INTEGER::siz
    
    ! Local
    DOUBLE PRECISION::div
    INTEGER::ii, jj, sIndx, eIndx
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize tMatx
    tMatx = 0.D0
    DO ii = 1, siz
        tMatx(ii, ii) = 1.D0
    END DO
    
    ! Perform factorization
    DO jj = 1, SIZE(convecIndex(:, 1))
        eIndx = convecIndex(jj, 1)
        sIndx = convecIndex(jj, 2)
        
        ! Now, from sIndx to eIndx - 1, eliminate the non diagonal values
        ! of row eIndx
        DO ii = sIndx, eIndx - 1
            ! First, divide this row by its diagonal value
            DO WHILE(ABS(1.D0 - aa(ii, ii)).GE.1.D-16)
                div = 1.D0/aa(ii, ii)
                aa(ii, ii) = aa(ii, ii)*div
                aa(eIndx, ii) = aa(eIndx, ii)*div
                tMatx(ii, ii) = tMatx(ii, ii)*div
            END DO
            
            ! Now eliminate
            aa(eIndx, eIndx) = aa(eIndx, eIndx) - aa(eIndx, ii)*aa(ii, eIndx)
            tMatx(ii, eIndx) = -tMatx(ii, ii)*aa(ii, eIndx)
            aa(ii, eIndx) = 0.D0
        END DO
        
        ! Make one the diagonal element of the convective shell
        DO WHILE(ABS(1.D0 - aa(eIndx, eIndx)).GE.1.D-16)
            div = 1.D0/aa(eIndx, eIndx)
            aa(eIndx, eIndx) = aa(eIndx, eIndx)*div
            tMatx(sIndx:eIndx, eIndx) = tMatx(sIndx:eIndx, eIndx)*div
        END DO
    END DO
    
END SUBROUTINE sparseOvMatrixGauss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine performs a backwards substitution over aa to solve xx    !!!
!!! having bb. After the subroutine exits, the solution will be in xx. The   !!!
!!! matrix tMatx is the transformation matrix for bb.                        !!!
!!! dimSiz is the additional dimension for bb (more than one bb at once).    !!!
!!! This subroutine is the complement to sparseGauss.                        !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -aa, the matrix in which perform the elimination.                        !!!
!!! -tMatx, the matrix holding all the transformations done during the       !!!
!!!          elimination. That way we can transform the independent vector.  !!!
!!! -bb, the array of independent terms arrays.                              !!!
!!! -xx, holds the solutions for the different bb on exit.                   !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the further index they affect.                     !!!
!!! -siz, aa size.                                                           !!!
!!! -dimSiz, the number of different bb we have to solve simultaneously.     !!!
!!!                                                                          !!!
!!! On output the solution is given in the array xx.                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dimBackSub(aa, tMatx, bb, xx, convecIndex, siz, dimSiz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::aa(:, :), tMatx(:, :), bb(:, :), xx(:, :)
    INTEGER::convecIndex(:, :), siz, dimSiz
    
    ! Local
    DOUBLE PRECISION::tot(dimSiz)
    INTEGER::ii, jj, jjindx, times
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! First perform bb transformation
    times = 1
    DO ii = 1, siz
        ! Identify if we are in a radiative or convective shell
        IF (ii.EQ.convecIndex(times, 1)) THEN
            tot = 0
            DO jj = convecIndex(times, 2), ii
                tot = tot + tMatx(jj, ii)*bb(:, jj)
            END DO
            
            xx(:, ii) = tot
            times = times + 1
        ELSE
            xx(:, ii) = tMatx(ii, ii)*bb(:, ii)
        END IF
    END DO
    
    ! Now solve for every convective shell
    times = 1
    DO WHILE (times.LE.SIZE(convecIndex(:, 1)))
        jjindx = convecIndex(times, 1)
        DO ii = jjindx - 1, convecIndex(times, 2), -1
            xx(:, ii) = xx(:, ii) - aa(jjindx, ii)*xx(:, jjindx)
        END DO
        
        times = times + 1
    END DO
    
END SUBROUTINE dimBackSub

END MODULE math_routines

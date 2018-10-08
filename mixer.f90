MODULE mixer
    USE MPI
    USE structures_mod
    USE math_routines
    IMPLICIT NONE
    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine mixes the convective zone and introduces it.             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -totShell, the size of intShell.                                         !!!
!!! -siz, the number of species (size of abundance dimension).               !!!
!!!                                                                          !!!
!!! On output intShell is filled and every convective zone is mixed          !!!
!!! instantaneously. The mixed values are in intShell and the common pointer !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE storeConvection(intShell, totShell, siz)
    IMPLICIT NONE
   
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:)
    INTEGER::totShell, siz
    
    ! Local
    TYPE(INTERSHELL), POINTER::commShell
    DOUBLE PRECISION::convdens(siz), convmass, dm
    INTEGER::ii
    LOGICAL::first
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize variables to shut up compiler
    NULLIFY(commShell)
    convdens = 0.D0; convmass = 0.D0
    
    ! Now do this for every convective region
    ii = 1
    DO
        first = .TRUE.
        DO WHILE (intShell(ii)%isConvective)
            IF (first) THEN
                first = .FALSE.
                
                ! Allocate the common shell
                ALLOCATE(commShell)
                ALLOCATE(commShell%dens(siz))
                
                ! Store starting mass and index
                commShell%ii0 = ii
                commShell%mass0 = intShell(ii)%mass0
                commShell%rad0 = intShell(ii)%rad0
                commShell%rho0 = intShell(ii)%rho0
                commShell%hp0 = intShell(ii)%hp0
                commShell%pres0 = intShell(ii)%pres0
                commShell%vel0 = intShell(ii)%vel0
                
                convdens = 0.D0; convmass = 0.D0
            END IF
            
            ! Put speeds
            IF (commShell%vel0.LT.1.D-50) commShell%vel0 = intShell(ii)%vel0
            IF (intShell(ii)%vel1.GT.1.D-50) commShell%vel1 = intShell(ii)%vel1
            
            ! Point with current shell
            intShell(ii)%convShell => commShell
            
            ! Store convdens and convmass
            dm = intShell(ii)%mass1 - intShell(ii)%mass0
            convdens = convdens + intShell(ii)%dens*dm
            convmass = convmass + dm
            
            ! Advance loop
            ii = ii + 1
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! If we had a convective region, introduce the last values
        ! and nullify
        IF (.NOT.first) THEN
            ! Correct index
            ii = ii - 1
            
            ! Introduce values
            commShell%ii1 = ii
            commShell%mass1 = intShell(ii)%mass1
            commShell%rad1 = intShell(ii)%rad1
            commShell%rho1 = intShell(ii)%rho1
            commShell%hp1 = intShell(ii)%hp1
            commShell%pres1 = intShell(ii)%pres1
            
            IF (intShell(ii)%vel1.GT.1.D-50) commShell%vel1 = intShell(ii)%vel1
            
            ! Beware of zero mass shells
            IF (convmass.GT.1.D-50) THEN
                commShell%dens = convdens/convmass
            ELSE
                commShell%dens = intShell(ii)%dens
            END IF
            
            ! Nullify for next convective region
            NULLIFY(commShell)
        END IF
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! One more loop to store the mixed abundance in the shells
    DO ii = 1, totShell - 1
        IF (intShell(ii)%isConvective) THEN
            intShell(ii)%dens = intShell(ii)%convShell%dens
        END IF
    END DO
    
END SUBROUTINE storeConvection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the advective overshooting limit.             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -p1indx, index for protons.                                              !!!
!!! -he4indx, index for alpha particles.                                     !!!
!!! -ovParam, the envelope overshooting parameter.                           !!!
!!! -ovPDCZParam, the PDCZ overshooting parameter.                           !!!
!!!                                                                          !!!
!!! On output firstOv holds the first overshooting index.                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE advOvLimit(intShell, firstOv, p1indx, he4indx, ovParam, ovPDCZParam)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), TARGET::intShell(:)
    DOUBLE PRECISION::ovParam, ovPDCZParam
    INTEGER::firstOv, p1indx, he4indx
    
    ! Local
    TYPE (INTERSHELL), POINTER::sh, ovShell(:)
    DOUBLE PRECISION::fthick, convecSpan, pb, vb, rb, dmb, ptild, omeg
    DOUBLE PRECISION::hps(2), ps(2), rs(2), rmed, kval, dms
    DOUBLE PRECISION::localOvParam
    INTEGER::ii, jj, jjndx, siz, totShell, nShells
    LOGICAL::inShell, entered
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get sizes
    siz = SIZE(intShell(1)%dens)
    totShell = SIZE(intShell) + 1
    
    ! Count the effective number of shells, treating convective regions as one
    ii = 1; nShells = 0
    DO
        ! Add one
        nShells = nShells + 1
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! Allocate ovShell
    ALLOCATE(ovShell(nShells))
    DO ii = 1, nShells
        ALLOCATE(ovShell(ii)%dens(siz))
    END DO
    
    ! Store values in ovShell
    ii = 1
    DO jj = 1, nShells
        ! Add this shell
        IF (intShell(ii)%isConvective) THEN
            sh => intShell(ii)%convShell
        ELSE
            sh => intShell(ii)
        END IF
        
        ! Copy every value
        ovShell(jj)%mass0 = sh%mass0
        ovShell(jj)%mass1 = sh%mass1
        ovShell(jj)%hp0 = sh%hp0
        ovShell(jj)%hp1 = sh%hp1
        ovShell(jj)%pres0 = sh%pres0
        ovShell(jj)%pres1 = sh%pres1
        ovShell(jj)%vel0 = sh%vel0
        ovShell(jj)%vel1 = sh%vel1
        ovShell(jj)%rad0 = sh%rad0
        ovShell(jj)%rad1 = sh%rad1
        ovShell(jj)%dens = sh%dens
        
        ! Copy convection flag
        ovShell(jj)%isConvective = intShell(ii)%isConvective
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
    END DO
    
    ! Initialize firstOv
    firstOv = totShell
    
    ! Go shell by shell checking
    localOvParam = ovParam
    DO ii = nShells, 1, -1
        ! Cycle if not convective shell
        IF (.NOT.ovShell(ii)%isConvective) CYCLE
        
        ! Check if we have to change the localOvParam, exit if in core
        IF (ovShell(ii)%dens(p1indx).LT.1.D-3) THEN
            IF (ovShell(ii)%dens(he4indx).LT.1.D-4) EXIT
            localOvParam = ovPDCZParam
        END IF
        
        ! Calculate and store convective zone mass
        dmb = ABS(ovShell(ii)%mass1 - ovShell(ii)%mass0)
        
        ! If dmb is too small ignore this convective zone
        ! This is done to avoid dividing by zero
        IF (dmb.LT.1.D-50) CYCLE
        
        ! Store convective (bubble) zone variables
        pb = ovShell(ii)%pres0
        vb = ovShell(ii)%vel0
        rb = ovShell(ii)%rad0
        
        ! Calculate convective zone span
        convecSpan = ovShell(ii)%rad1 - ovShell(ii)%rad0
        
        ! Set flag to true
        inShell = .TRUE.
        
        DO jj = 1, nShells
            ! Calculate the index
            jjndx = ii - jj
            
            ! Borders and convective regions
            IF (((jjndx).LT.1).OR.(jjndx.GT.nShells)) THEN
                inShell = .FALSE.
            ELSE IF (ovShell(jjndx)%isConvective) THEN
                inShell = .FALSE.
            END IF
            
            ! Exit if no more overshooting is to be applied
            IF (.NOT.inShell) THEN
                firstOv = jjndx
                EXIT
            END IF
            
            ! Calculate shell mass and cycle if zero mass shell
            dms = ABS(ovShell(jjndx)%mass1 - ovShell(jjndx)%mass0)
            IF (dms.LT.1.D-50) CYCLE
            
            ! Store shell variables
            hps(1) = ovShell(jjndx)%hp0
            hps(2) = ovShell(jjndx)%hp1
            
            ps(1) = ovShell(jjndx)%pres0
            ps(2) = ovShell(jjndx)%pres1
            
            rs(1) = ovShell(jjndx)%rad0
            rs(2) = ovShell(jjndx)%rad1
            
            ! Calculate the speed
            rmed = SUM(rs)*0.5D0
            ptild = SUM(ps)/pb*0.5D0
            
            fthick = convecSpan/(SUM(hps)*0.5D0)
            IF (fthick.GT.1.D0) fthick = 1.D0
            
            omeg = 1.D0/(localOvParam*fthick)
            IF (ptild.GT.1.D0) omeg = -omeg
            
            kval = vb*(ptild**omeg)
            
            ! If the speed (vb*ptild**omeg) < 1.d-10 that is our limit
            IF (kval.LT.1.D-10) THEN
                inShell = .FALSE.
                firstOv = jjndx
                EXIT
            END IF
        END DO
    END DO
    
    ! Keep borders in bounds
    IF (firstOv.LT.1) firstOv = 1
    IF (firstOv.GT.nShells) firstOv = nShells
    
    ! Deallocate ovShell
    DO ii = 1, nShells
        DEALLOCATE(ovShell(ii)%dens)
    END DO
    DEALLOCATE(ovShell)
    
END SUBROUTINE advOvLimit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine calculates the diffusive overshooting limit.             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -p1indx, index for protons.                                              !!!
!!! -he4indx, index for alpha particles.                                     !!!
!!! -ovParam, the envelope overshooting parameter.                           !!!
!!! -ovPDCZParam, the PDCZ overshooting parameter.                           !!!
!!!                                                                          !!!
!!! On output firstOv holds the first overshooting index.                    !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE diffOvLimit(intShell, firstOv, p1indx, he4indx, ovParam, ovPDCZParam)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), TARGET::intShell(:)
    DOUBLE PRECISION::ovParam, ovPDCZParam
    INTEGER::firstOv, p1indx, he4indx
    
    ! Local
    TYPE (INTERSHELL), POINTER::sh, ovShell(:)
    DOUBLE PRECISION::fthick, convecSpan, pb, vb, rb, dmb, ptild, omeg, hps, ps
    DOUBLE PRECISION::rs, kval, dms, rhob, rhos, hpb, localOvParam
    DOUBLE PRECISION, PARAMETER::pi = 4*ATAN(1.D0), msun = 1.9984D30
    INTEGER::ii, jj, jjndx, siz, totShell, nShells
    LOGICAL::inShell, entered
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get sizes
    siz = SIZE(intShell(1)%dens)
    totShell = SIZE(intShell) + 1
    
    ! Count the effective number of shells, treating convective regions as one
    ii = 1; nShells = 0
    DO
        ! Add one
        nShells = nShells + 1
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! Allocate ovShell
    ALLOCATE(ovShell(nShells))
    DO ii = 1, nShells
        ALLOCATE(ovShell(ii)%dens(siz))
    END DO
    
    ! Store values in ovShell
    ii = 1
    DO jj = 1, nShells
        ! Add this shell
        IF (intShell(ii)%isConvective) THEN
            sh => intShell(ii)%convShell
        ELSE
            sh => intShell(ii)
        END IF
        
        ! Copy every value
        ovShell(jj)%mass0 = sh%mass0
        ovShell(jj)%mass1 = sh%mass1
        ovShell(jj)%rho0 = sh%rho0
        ovShell(jj)%rho1 = sh%rho1
        ovShell(jj)%hp0 = sh%hp0
        ovShell(jj)%hp1 = sh%hp1
        ovShell(jj)%pres0 = sh%pres0
        ovShell(jj)%pres1 = sh%pres1
        ovShell(jj)%vel0 = sh%vel0
        ovShell(jj)%vel1 = sh%vel1
        ovShell(jj)%rad0 = sh%rad0
        ovShell(jj)%rad1 = sh%rad1
        ovShell(jj)%dens = sh%dens
        
        ! Copy convection flag
        ovShell(jj)%isConvective = intShell(ii)%isConvective
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
    END DO
    
    ! Initialize firstOv
    firstOv = totShell
    
    ! Go shell by shell checking
    localOvParam = ovParam
    DO ii = nShells, 1, -1
        ! Cycle if not convective shell
        IF (.NOT.ovShell(ii)%isConvective) CYCLE
        
        ! Check if we have to change the localOvParam, exit if in core
        IF (ovShell(ii)%dens(p1indx).LT.1.D-3) THEN
            IF (ovShell(ii)%dens(he4indx).LT.1.D-4) EXIT
            localOvParam = ovPDCZParam
        END IF
        
        ! Calculate and store convective zone mass
        dmb = ABS(ovShell(ii)%mass1 - ovShell(ii)%mass0)
        
        ! If dmb is too small ignore this convective zone
        ! This is done to avoid dividing by zero
        IF (dmb.LT.1.D-50) CYCLE
        
        ! Store convective (bubble) zone variables
        pb = ovShell(ii)%pres0
        rhob = ovShell(ii)%rho0/msun ! In solar masses
        vb = ovShell(ii)%vel0
        hpb = ovShell(ii)%hp0
        rb = ovShell(ii)%rad0
        
        ! Calculate convective zone span
        convecSpan = ovShell(ii)%rad1 - rb
        
        ! Set flag to true
        inShell = .TRUE.
        
        DO jj = 1, nShells
            ! Calculate the index
            jjndx = ii - jj
            
            ! Borders and convective regions
            IF (((jjndx).LT.1).OR.(jjndx.GT.nShells)) THEN
                inShell = .FALSE.
            ELSE IF (ovShell(jjndx)%isConvective) THEN
                inShell = .FALSE.
            END IF
            
            ! Exit if no more overshooting is to be applied
            IF (.NOT.inShell) THEN
                firstOv = jjndx
                EXIT
            END IF
            
            ! Calculate shell mass and cycle if zero mass shell
            dms = ABS(ovShell(jjndx)%mass1 - ovShell(jjndx)%mass0)
            IF (dms.LT.1.D-50) CYCLE
            
            ! Store shell variables
            hps = ovShell(jjndx)%hp0
            ps = ovShell(jjndx)%pres0
            rhos = ovShell(jjndx)%rho0/msun ! In solar masses
            rs = ovShell(jjndx)%rad0
            
            ! Calculate the adimensional pressure
            ptild = ps/pb
            
            fthick = convecSpan/hps
            IF (fthick.GT.1.D0) fthick = 1.D0
            
            omeg = 1.D0/(localOvParam*fthick)
            IF (ptild.GT.1.D0) omeg = -omeg
            
            kval = (4*pi*rs*rs*rhos)**2
            kval = kval*vb*(ptild**omeg)*0.1*hpb/3.D0
            
            ! If the diffusion coefficient is lower than 10^-30
            ! which would mean a tau of 10^16 seconds for a
            ! characteristic lenght of 10^-7, we consider it done
            IF (kval.LT.1.D-30) THEN
                inShell = .FALSE.
                firstOv = jjndx
                EXIT
            END IF
        END DO
    END DO
    
    ! Keep borders in bounds
    IF (firstOv.LT.1) firstOv = 1
    IF (firstOv.GT.nShells) firstOv = nShells
    
    ! Deallocate ovShell
    DO ii = 1, nShells
        DEALLOCATE(ovShell(ii)%dens)
    END DO
    DEALLOCATE(ovShell)
    
END SUBROUTINE diffOvLimit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine does a convective mixing from values in intShell.        !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -totShell, the size of intShell.                                         !!!
!!! -siz, the number of species (size of abundance dimension).               !!!
!!!                                                                          !!!
!!! On output all convective zones are instantanously mixed. The values are  !!!
!!! in the common pointer of convective shells.                              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE mixConvection(intShell, totShell, siz)
    IMPLICIT NONE
   
    ! Input
    TYPE (INTERSHELL)::intShell(:)
    INTEGER::totShell, siz
    
    ! Local
    TYPE(INTERSHELL), POINTER::commShell
    DOUBLE PRECISION::convdens(siz), convmass, dm
    INTEGER::ii
    LOGICAL::first
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Initialize
    NULLIFY(commShell)
    
    ! Do this for every convective region
    ii = 1
    DO
        first = .TRUE.
        DO WHILE (intShell(ii)%isConvective)
            IF (first) THEN
                first = .FALSE.
                convdens = 0.D0; convmass = 0.D0
                commShell => intShell(ii)%convShell
            END IF
            
            ! Store convdens and convmass
            dm = intShell(ii)%mass1 - intShell(ii)%mass0
            convdens = convdens + intShell(ii)%dens*dm
            convmass = convmass + dm
            
            ! Advance loop
            ii = ii + 1
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index and introduce values
        IF (.NOT.first) THEN
            ii = ii - 1
            
            ! Beware of zero mass shells
            IF (convmass.GT.1.D-50) THEN
                commShell%dens = convdens/convmass
            ELSE
                commShell%dens = intShell(ii)%dens
            END IF
        END IF
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! Return the mixed values to the intShell
    DO ii = 1, totShell - 1
        IF (intShell(ii)%isConvective) THEN
            intShell(ii)%dens = intShell(ii)%convShell%dens
        END IF
    END DO
    
END SUBROUTINE mixConvection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine does the final mixing from values in the intShell and    !!!
!!! cleans the convective shell.                                             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -totShell, the size of intShell.                                         !!!
!!! -siz, the number of species (size of abundance dimension).               !!!
!!!                                                                          !!!
!!! On output all convective zones are instantanously mixed and all pointers !!!
!!! are nullified. The mixed values are in intShell.                         !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE cleanConvection(intShell, totShell, siz)
    IMPLICIT NONE
   
    ! Input
    TYPE (INTERSHELL), POINTER::intShell(:)
    INTEGER::totShell, siz
    
    ! Local
    INTEGER::ii, iiend
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Start by mixing the convective region
    CALL mixConvection(intShell, totShell, siz)
    
    ! Store the mixed abundance in the shells
    DO ii = 1, totShell - 1
        IF (intShell(ii)%isConvective) THEN
            intShell(ii)%dens = intShell(ii)%convShell%dens
        END IF
    END DO
    
    ! Now deallocate and nullify the convShell common pointer
    ii = 1
    DO ii = 1, totShell - 1
        IF (intShell(ii)%isConvective) THEN
            iiend = intShell(ii)%convShell%ii1
            
            ! Deallocate
            IF (ii.EQ.iiend) THEN
                DEALLOCATE(intShell(ii)%convShell%dens)
                DEALLOCATE(intShell(ii)%convShell)
            END IF
            
            ! Nullify
            NULLIFY(intShell(ii)%convShell)
        END IF
    END DO
    
END SUBROUTINE cleanConvection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the overshooting matrix.                         !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -ovMatrix, the overshooting matrix.                                      !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -ovShell, the shell arrays for the overshooting matrix.                  !!!
!!! -nShells, size of ovShell.                                               !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -p1indx, index for protons.                                              !!!
!!! -he4indx, index for alpha particles.                                     !!!
!!! -ovParam, the envelope overshooting parameter.                           !!!
!!! -ovPDCZParam, the PDCZ overshooting parameter.                           !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the furthermost index they affect.                 !!!
!!! -performOv, bool to know if we should perform overshooting or not.       !!!
!!!                                                                          !!!
!!! On output all extra-mixed values are stored in intShell, even the        !!!
!!! convective ones. The common convective pointer is outdated.              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createOvMatrix(ovMatrix, intShell, ovShell, nShells, firstOv, &
                          p1indx, he4indx, ovParam, ovPDCZParam, convecIndex, &
                          performOv)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::ovShell(:)
    TYPE (INTERSHELL), TARGET::intShell(:)
    DOUBLE PRECISION, POINTER::ovMatrix(:, :)
    DOUBLE PRECISION::ovParam, ovPDCZParam
    INTEGER::nShells, firstOv, p1indx, he4indx
    INTEGER, POINTER::convecIndex(:, :)
    LOGICAL::performOv
    
    ! Local
    TYPE (INTERSHELL), POINTER::sh
    DOUBLE PRECISION, POINTER::holdOvMatrix(:, :)
    DOUBLE PRECISION::fthick, convecSpan, pb(2), vb(2), rb(2), dmb, ptild, omeg
    DOUBLE PRECISION::hps(2), ps(2), rs(2), rmed, tauval, kval, dms, timeVal
    DOUBLE PRECISION::preRm, localOvParam
    INTEGER::ii, jj, mm, jjndx, sig(2), siz, totShell, lastOv, times, sIndx
    INTEGER::eIndx, locFirstOv
    LOGICAL::inShell, entered
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get sizes
    siz = SIZE(intShell(1)%dens)
    totShell = SIZE(intShell) + 1
    
    ! Count the effective number of shells, treating convective regions as one
    ii = 1; nShells = 0
    DO
        ! Add one
        nShells = nShells + 1
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! Allocate ovShell and ovMatrix
    ALLOCATE(ovShell(nShells), ovMatrix(nShells, nShells))
    DO ii = 1, nShells
        ALLOCATE(ovShell(ii)%dens(siz))
    END DO
    
    ! Store values in ovShell
    ii = 1
    DO jj = 1, nShells
        ! Add this shell
        IF (intShell(ii)%isConvective) THEN
            sh => intShell(ii)%convShell
        ELSE
            sh => intShell(ii)
        END IF
        
        ! Copy every value
        ovShell(jj)%mass0 = sh%mass0
        ovShell(jj)%mass1 = sh%mass1
        ovShell(jj)%hp0 = sh%hp0
        ovShell(jj)%hp1 = sh%hp1
        ovShell(jj)%pres0 = sh%pres0
        ovShell(jj)%pres1 = sh%pres1
        ovShell(jj)%vel0 = sh%vel0
        ovShell(jj)%vel1 = sh%vel1
        ovShell(jj)%rad0 = sh%rad0
        ovShell(jj)%rad1 = sh%rad1
        ovShell(jj)%dens = sh%dens
        
        ! Copy convection flag
        ovShell(jj)%isConvective = intShell(ii)%isConvective
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
    END DO
    
    ! This matrix will hold the coefficients for the overshooting transformation
    ! These coefficients allow us to write the system
    ! d(X_s)/dt = Ksc*(X_c - X_s)
    ! d(X_c)/dt = -ms/mc*Ksc(X_c - X_s)
    
    ovMatrix = 0.D0
    
    ! Sign array
    sig = (/-1, 1/)
    
    ! Go shell by shell creating the coefficients
    localOvParam = ovParam
    DO ii = nShells, 1, -1
        ! Cycle if not convective shell
        IF (.NOT.ovShell(ii)%isConvective) CYCLE
        
        ! Check if we have to change the localOvParam, exit if in core
        IF (ovShell(ii)%dens(p1indx).LT.1.D-3) THEN
            IF (ovShell(ii)%dens(he4indx).LT.1.D-4) EXIT
            localOvParam = ovPDCZParam
        END IF
        
        ! Calculate and store convective zone mass
        dmb = ABS(ovShell(ii)%mass1 - ovShell(ii)%mass0)
        
        ! If dmb is too small ignore this convective zone
        ! This is done to avoid dividing by zero
        IF (dmb.LT.1.D-50) CYCLE
        
        ! Store convective (bubble) zone variables
        pb(1) = ovShell(ii)%pres0
        pb(2) = ovShell(ii)%pres1
        vb(1) = ovShell(ii)%vel0
        vb(2) = ovShell(ii)%vel1
        rb(1) = ovShell(ii)%rad0
        rb(2) = ovShell(ii)%rad1
        
        ! Calculate convective zone span
        convecSpan = rb(2) - rb(1)
        
        ! Apply overshooting
        DO mm = 1, 2
            ! Apply it only downwards
            IF (mm.EQ.2) CYCLE ! WARNING:
                               ! The code can only go downwards for now
            
            ! Set flag to true
            inShell = .TRUE.
            
            ! Set to 0 the timeVal and to preRm rb
            timeVal = 0.D0
            preRm = rb(mm)
            
            DO jj = 1, nShells
                ! Calculate the index
                jjndx = ii + sig(mm)*jj
                
                ! Borders and convective regions
                IF (((jjndx).LT.1).OR.(jjndx.GT.nShells)) THEN
                    inShell = .FALSE.
                ELSE IF (ovShell(jjndx)%isConvective) THEN
                    inShell = .FALSE.
                END IF
                
                ! Exit if no more overshooting is to be applied
                IF (.NOT.inShell) EXIT
                
                ! Calculate shell mass and cycle if zero mass shell
                dms = ABS(ovShell(jjndx)%mass1 - ovShell(jjndx)%mass0)
                IF (dms.LT.1.D-50) CYCLE
                
                ! Store shell variables
                hps(1) = ovShell(jjndx)%hp0
                hps(2) = ovShell(jjndx)%hp1
                
                ps(1) = ovShell(jjndx)%pres0
                ps(2) = ovShell(jjndx)%pres1
                
                rs(1) = ovShell(jjndx)%rad0
                rs(2) = ovShell(jjndx)%rad1
                
                ! Calculate the speed
                rmed = SUM(rs)*0.5D0
                ptild = SUM(ps)/pb(mm)*0.5D0
                
                fthick = convecSpan/(SUM(hps)*0.5D0)
                IF (fthick.GT.1.D0) fthick = 1.D0
                
                omeg = 1.D0/(localOvParam*fthick)
                IF (ptild.GT.1.D0) omeg = -omeg
                
                kval = vb(mm)*(ptild**omeg)
                
                ! If the speed (vb*ptild**omeg) is lesser than
                ! a certain value in cm/s, make it zero
                IF (kval.LT.1.D-10) THEN
                    inShell = .FALSE.
                    kval = 0.D0
                    EXIT
                ELSE
                    ! Calculate kval and introduce it
                    timeVal = timeVal + ABS(rmed - preRm)/kval
                    kval = 1/timeVal
                    tauval = dms/dmb*kval
                    
                    ! Introduce in shell row
                    ovMatrix(jjndx, ii) = ovMatrix(jjndx, ii) + kval
                    ovMatrix(jjndx, jjndx) = ovMatrix(jjndx, jjndx) - kval
                    
                    ! Introduce in convective zone row
                    ovMatrix(ii, ii) = ovMatrix(ii, ii) - tauval
                    ovMatrix(ii, jjndx) = ovMatrix(ii, jjndx) + tauval
                    
                    ! Update preRm
                    preRm = rmed
                END IF
            END DO
        END DO
    END DO
    
    ! Check matrix to know if perform overshooting
    performOv = .FALSE.
    DO ii = 1, nShells
        IF (ABS(ovMatrix(ii, ii)).GT.1.D-300) THEN
            performOv = .TRUE.
            EXIT
        END IF
    END DO
    
    ! If perform overshooting, reduce ovMatrix
    locFirstOv = 0; lastOv = 0
    IF (performOv) THEN
        ! Get first and last ov shells
        DO ii = 1, nShells
            IF (ABS(ovMatrix(ii, ii)).GT.1.D-300) THEN
                IF (locFirstOv.EQ.0) THEN
                    locFirstOv = ii
                ELSE
                    lastOv = ii
                END IF
            END IF
        END DO
        
        IF (locFirstOv.LT.firstOv) THEN
            PRINT*, "firstOv > locFirstOv, stopping"
            PRINT*, firstOv, locFirstOv
            STOP
        END IF
        
        ! Store big matrix
        holdOvMatrix => ovMatrix
        NULLIFY(ovMatrix)
        
        ! Redefine nShells
        nShells = lastOv - firstOv + 1
        
        ! Allocate and initialize ovMatrix
        ALLOCATE(ovMatrix(nShells, nShells))
        ovMatrix = 0
        
        ! Now copy the relevant rows of holdOvMatrix
        DO ii = 1, nShells
            ovMatrix(ii, :) = holdOvMatrix(ii - 1 + firstOv, firstOv:lastOv)
        END DO
        
        ! Deallocate old matrix
        DEALLOCATE(holdOvMatrix)
        NULLIFY(holdOvMatrix)
        
        ! Create convecIndex array
        times = 0
        DO ii = 2, nShells
            IF (ABS(ovMatrix(ii - 1, ii)).GT.1.D-300) times = times + 1
        END DO
        
        ! Allocate
        ALLOCATE(convecIndex(times, 2))
        
        sIndx = 1; eIndx = 1; times = 1
        DO WHILE (sIndx.LT.nShells)
            ! Look for eIndx. The eIndx will be the location
            ! of the next convective shell
            DO ii = sIndx + 1, nShells
                IF (ABS(ovMatrix(ii - 1, ii)).GT.1.D-300) THEN
                    eIndx = ii
                    EXIT
                END IF
            END DO
            
            ! Put information into convecIndex
            convecIndex(times, 1) = eIndx
            convecIndex(times, 2) = sIndx
            times = times + 1
            
            ! Now look for next sIndx, where the matrix
            ! diagonal is greater than zero
            sIndx = eIndx + 1
            DO ii = sIndx, nShells
                IF (ABS(ovMatrix(ii, ii)).GT.1.D-300) THEN
                    sIndx = ii
                    EXIT
                END IF
            END DO
        END DO
    END IF
    
END SUBROUTINE createOvMatrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine creates the overshooting matrix.                         !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -ovMatrix, a matrix holding the overshooting arrays (a, b and c).        !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -ovShell, the shell arrays for the overshooting matrix.                  !!!
!!! -nShells, size of ovShell.                                               !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -p1indx, index for protons.                                              !!!
!!! -he4indx, index for alpha particles.                                     !!!
!!! -ovParam, the envelope overshooting parameter.                           !!!
!!! -ovPDCZParam, the PDCZ overshooting parameter.                           !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the furthermost index they affect.                 !!!
!!! -dx, delta x, with dx(j) defined as x_{j + 1} - x_j                      !!!
!!! -performOv, bool to know if we should perform overshooting or not.       !!!
!!!                                                                          !!!
!!! On output all extra-mixed values are stored in intShell, even the        !!!
!!! convective ones. The common convective pointer is outdated.              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE createOvArrays(ovMatrix, intShell, ovShell, nShells, firstOv, &
                          p1indx, he4indx, ovParam, ovPDCZParam, convecIndex, &
                          dx, performOv)
    IMPLICIT NONE
    
    ! Input
    TYPE (INTERSHELL), POINTER::ovShell(:)
    TYPE (INTERSHELL), TARGET::intShell(:)
    DOUBLE PRECISION, POINTER::ovMatrix(:, :), dx(:)
    DOUBLE PRECISION::ovParam, ovPDCZParam
    INTEGER::nShells, firstOv, p1indx, he4indx
    INTEGER, POINTER::convecIndex(:, :)
    LOGICAL::performOv
    
    ! Local
    TYPE (INTERSHELL), POINTER::sh
    DOUBLE PRECISION, POINTER::holdOvMatrix(:, :), holdDx(:)
    DOUBLE PRECISION::fthick, convecSpan, pb, rb, rhob, vb, dmb, ptild, omeg
    DOUBLE PRECISION::hps, ps, rs, rhos, kval, dms, hpb, localOvParam, radRhob
    DOUBLE PRECISION, PARAMETER::pi = 4*ATAN(1.D0), msun = 1.9984D30
    INTEGER, POINTER::holdConvecIndex(:, :)
    INTEGER::ii, jj, mm, jjndx, sig(2), siz, totShell, times
    LOGICAL::inShell, entered
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get sizes
    siz = SIZE(intShell(1)%dens)
    totShell = SIZE(intShell) + 1
    
    ! Count the effective number of shells, treating convective regions as one
    ii = 1; nShells = 0
    DO
        ! Add one
        nShells = nShells + 1
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
        IF (ii.GE.totShell) EXIT
    END DO
    
    ! Allocate ovShell and ovMatrix
    ALLOCATE(ovShell(nShells), ovMatrix(nShells, 4))
    ALLOCATE(convecIndex(nShells, 2), dx(nShells))
    DO ii = 1, nShells
        ALLOCATE(ovShell(ii)%dens(siz))
    END DO
    
    ! Store values in ovShell
    ii = 1; dx = 1.D0
    DO jj = 1, nShells
        ! Add this shell
        IF (intShell(ii)%isConvective) THEN
            sh => intShell(ii)%convShell
        ELSE
            sh => intShell(ii)
        END IF
        
        ! Copy every value
        ovShell(jj)%mass0 = sh%mass0
        ovShell(jj)%mass1 = sh%mass1
        ovShell(jj)%rho0 = sh%rho0
        ovShell(jj)%rho1 = sh%rho1
        ovShell(jj)%hp0 = sh%hp0
        ovShell(jj)%hp1 = sh%hp1
        ovShell(jj)%pres0 = sh%pres0
        ovShell(jj)%pres1 = sh%pres1
        ovShell(jj)%vel0 = sh%vel0
        ovShell(jj)%vel1 = sh%vel1
        ovShell(jj)%rad0 = sh%rad0
        ovShell(jj)%rad1 = sh%rad1
        ovShell(jj)%dens = sh%dens
        
        ! Copy convection flag
        ovShell(jj)%isConvective = intShell(ii)%isConvective
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
    END DO
    
    ! Create dx
    DO ii = 1, nShells
        dx(ii) = ovShell(ii)%mass1 - ovShell(ii)%mass0
    END DO
    
    ! The ovMatrix holds the coefficients a, b and c for the solution of the
    ! diffusive overshooting with the backwards euler method
    
    ! ovMatrix(:, 1) = a(:); ovMatrix(:, 2) = b(:); ovMatrix(:, 3) = c and
    ! ovMatrix(:, 4) = k(:) <-- diffusion coefficient.
    ovMatrix = 0.D0
    
    ! Sign array
    sig = (/-1, 1/)
    
    ! Go shell by shell creating the coefficients
    times = 0; convecIndex = 0
    localOvParam = ovParam
    DO ii = nShells, 1, -1
        
        ! Cycle if not convective shell
        IF (.NOT.ovShell(ii)%isConvective) CYCLE
        
        ! Check if we have to change the localOvParam, exit if in core
        IF (ovShell(ii)%dens(p1indx).LT.1.D-3) THEN
            IF (ovShell(ii)%dens(he4indx).LT.1.D-4) EXIT
            localOvParam = ovPDCZParam
        END IF
        
        ! Count this as a convective zone
        times = times + 1
        
        ! Calculate and store convective zone mass
        dmb = ABS(ovShell(ii)%mass1 - ovShell(ii)%mass0)
        
        ! If dmb is too small ignore this convective zone
        ! This is done to avoid dividing by zero
        IF (dmb.LT.1.D-50) CYCLE
        
        ! Store convective (bubble) zone variables
        pb = ovShell(ii)%pres0
        rhob = ovShell(ii)%rho0/msun ! In solar masses
        vb = ovShell(ii)%vel0
        hpb = ovShell(ii)%hp0
        rb = ovShell(ii)%rad0
        
        ! Calculate convective zone span
        convecSpan = ovShell(ii)%rad1 - rb
        
        ! Get diffusive coefficient
        radRhob = (4*pi*rb*rb*rhob)**2
        radRhob = radRhob*0.1*hpb/3.D0
        
        ! Add it to the matrix. We do not care about the coefficient at ii + 1
        ! because the flux 0 condition takes care of it.
        ovMatrix(ii, 4) = vb*radRhob
        
        ! Fill convecIndex
        convecIndex(times, 1) = ii
        
        ! Apply overshooting
        DO mm = 1, 2
            ! Apply it only downwards
            IF (mm.EQ.2) CYCLE ! WARNING:
                               ! The code can only go downwards for now
            
            ! Set flag to true
            inShell = .TRUE.
            
            DO jj = 1, nShells
                ! Calculate the index
                jjndx = ii + sig(mm)*jj
                
                ! Borders and convective regions
                IF (((jjndx).LT.1).OR.(jjndx.GT.nShells)) THEN
                    inShell = .FALSE.
                ELSE IF (ovShell(jjndx)%isConvective) THEN
                    inShell = .FALSE.
                END IF
                
                ! Exit if no more overshooting is to be applied
                IF (.NOT.inShell) EXIT
                
                ! Calculate shell mass and cycle if zero mass shell
                dms = ABS(ovShell(jjndx)%mass1 - ovShell(jjndx)%mass0)
                IF (dms.LT.1.D-50) CYCLE
                
                ! Store shell variables
                hps = ovShell(jjndx)%hp0
                ps = ovShell(jjndx)%pres0
                rhos = ovShell(jjndx)%rho0/msun ! In solar masses
                rs = ovShell(jjndx)%rad0
                
                ! Calculate the adimensional pressure
                ptild = ps/pb
                
                fthick = convecSpan/hps
                IF (fthick.GT.1.D0) fthick = 1.D0
                
                omeg = 1.D0/(localOvParam*fthick)
                IF (ptild.GT.1.D0) omeg = -omeg
                
                kval = (4*pi*rs*rs*rhos)**2
                kval = kval*vb*(ptild**omeg)*0.1*hpb/3.D0
                
                ! If the diffusion coefficient is lower than 10^-30
                ! which would mean a tau of 10^16 seconds for a
                ! characteristic lenght of 10^-7, we consider it done
                IF (kval.LT.1.D-30) THEN
                    inShell = .FALSE.
                    kval = 0.D0
                    
                    ! Check that there is at least one shell
                    ! Otherwise, eliminate this convective shell from the
                    ! calculations
                    IF (jj.EQ.1) THEN
                        convecIndex(times, :) = 0
                        times = times - 1
                    END IF
                    
                    EXIT
                ELSE
                    ! Add the diffusion coefficient to the array
                    ovMatrix(jjndx, 4) = kval
                    
                    ! Put this shell into convecIndex
                    convecIndex(times, 2) = jjndx
                END IF
            END DO
        END DO
    END DO
    
    ! Check if we should perform overshooting
    IF (times.GT.0) THEN
        performOv = .TRUE.
    ELSE
        performOv = .FALSE.
    END IF
    
    ! If performing overshooting, shrink arrays
    IF (performOv) THEN
        holdOvMatrix => ovMatrix
        holdConvecIndex => convecIndex
        holdDx => dx
        
        ! Take first and last indices
        ii = holdConvecIndex(times, 2)
        jj = holdConvecIndex(1, 1)
        
        IF (ii.LT.firstOv) THEN
            PRINT*, "firstOv > locFirstOv, stopping"
            PRINT*, firstOv, ii
            STOP
        END IF
        
        ! Now allocate new ovMatrix and convecIndex
        ! Redefine nShells first
        nShells = jj - firstOv + 1
        ALLOCATE(convecIndex(times, 2), ovMatrix(nShells, 4), dx(nShells))
        
        ovMatrix = holdOvMatrix(firstOv:jj, :)
        dx = holdDx(firstOv:jj)
        convecIndex = holdConvecIndex(1:times, :) - firstOv + 1
        
        DEALLOCATE(holdOvMatrix, holdConvecIndex)
        NULLIFY(holdOvMatrix, holdConvecIndex)
    END IF
    
    ! If perform overshooting, create coefficients
    
    ! The schematic for the coefficients according to the diffusion
    ! value k is (python notation):
    !
    ! a[1:-1] = k[2:]/(dx[1:-1]*dxt[2:])
    ! c[1:-1] = k[1:-1]/(dx[1:-1]*dxt[1:-1])
    ! b = a + c
    
    IF (performOv) THEN
        ! Create "a" coefficient
        DO ii = 1, nShells - 1
            ovMatrix(ii, 1) = ovMatrix(ii + 1, 4)/(dx(ii)*SUM(dx(ii:ii+1)))
            ovMatrix(ii, 1) = ovMatrix(ii, 1)*2
        END DO
        
        ! Now "c" coefficient
        DO ii = 2, nShells
            ovMatrix(ii, 3) = ovMatrix(ii, 4)/(dx(ii)*SUM(dx(ii-1:ii)))
            ovMatrix(ii, 3) = ovMatrix(ii, 3)*2
        END DO
        
        ! Finally "b" coefficient
        ovMatrix(:, 2) = ovMatrix(:, 1) + ovMatrix(:, 3)
    END IF
    
END SUBROUTINE createOvArrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine applies overshooting from values in intShell.            !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -intShell, the inter-shell array.                                        !!!
!!! -totShell, the size of intShell.                                         !!!
!!! -ovMatrix, the overshooting matrix.                                      !!!
!!! -ovShell, the shell arrays for the overshooting matrix.                  !!!
!!! -nShells, reduced size of ovMatrix.                                      !!!
!!! -firstOv, first shell affected by overshooting.                          !!!
!!! -dt, overshooting timestep.                                              !!!
!!! -eps, relative accuracy.                                                 !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the furthermost index they affect.                 !!!
!!! -yscale, threshold for zero abundances.                                  !!!
!!! -siz, the number of species (size of abundance dimension).               !!!
!!! -ovMode, the overshooting mode (advective or diffusive).                 !!!
!!! -nProc, total number of threads.                                         !!!
!!! -rank, this thread index.                                                !!!
!!!                                                                          !!!
!!! On output all extra-mixed values are stored in intShell, even the        !!!
!!! convective ones. The common convective pointer is outdated.              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE applyOvershooting(intShell, totShell, ovMatrix, ovShell, nShells, &
                            firstOv, dt, eps, convecIndex, yscale, siz, &
                            ovMode, nProc, rank)
    IMPLICIT NONE
    
    ! MPI variables
    INTEGER::ierror
    
    ! Input
    TYPE (INTERSHELL), TARGET::intShell(:), ovShell(:)
    DOUBLE PRECISION::ovMatrix(:, :), dt, eps, yscale
    INTEGER::totShell, nShells, firstOv, convecIndex(:, :), siz, nProc, rank
    CHARACTER(20)::ovMode
    
    ! Local
    TYPE (INTERSHELL), POINTER::sh
    DOUBLE PRECISION, ALLOCATABLE::redSol(:, :)
    DOUBLE PRECISION::prevSol(nShells, siz), dtot
    INTEGER::ii, jj, kk, redSiz, startIndx, endIndx, calcProc
    LOGICAL::entered
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Calculate the reduced size
    IF (nProc.GT.1) THEN
        IF (nProc.GE.siz) THEN
            redSiz = 1
            IF (rank.GE.siz) redSiz = 0
        ELSE
            redSiz = siz/nProc
            
            IF (rank.EQ.(nProc - 1)) THEN
                redSiz = siz - redSiz*(nProc - 1)
            END IF
        END IF
        
        ! Start and end indices
        IF (rank.NE.(nProc - 1)) THEN
            startIndx = rank*redSiz + 1
            endIndx = (rank + 1)*redSiz
        ELSE
            startIndx = siz - redSiz + 1
            endIndx = siz
        END IF
    ELSE
        startIndx = 1
        endIndx = siz
        redSiz = siz
    END IF
    
    ! Allocate prevSol
    ALLOCATE(redSol(nShells, redSiz))
    
    ! Add updated values to ovShell
    ii = 1
    DO jj = 1, SIZE(ovShell)
        ! Add this shell
        IF (intShell(ii)%isConvective) THEN
            sh => intShell(ii)%convShell
        ELSE
            sh => intShell(ii)
        END IF
        
        ! Copy abundance (only thing that has changed)
        ovShell(jj)%dens = sh%dens
        
        ! Skip convectives
        entered = .FALSE.
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            ii = ii + 1
            
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index
        IF (entered) ii = ii - 1
        
        ii = ii + 1
    END DO
    
    ! Intialize prevSol and redSol arrays
    DO jj = 1, nShells
        prevSol(jj, :) = ovShell(firstOv + jj - 1)%dens
        redSol(jj, :) = prevSol(jj, startIndx:endIndx)
    END DO
    
    ! Solve for dtot
    dtot = dt
    IF (ovMode.EQ."advective") THEN
        CALL solveAdvOv(ovMatrix, redSol, nShells, dtot, convecIndex, yscale, &
                        eps, redSiz)
    ELSE IF (ovMode.EQ."diffusive") THEN
        CALL solveDifOv(ovMatrix, redSol, nShells, dtot, yscale, eps, redSiz)
    END IF
    
    ! Copy values to prevSol
    DO jj = 1, nShells
        prevSol(jj, startIndx:endIndx) = redSol(jj, :)
    END DO
    
    ! Broadcast values
    IF (nProc.GT.1) THEN
        DO jj = 1, siz
            ! Calculate broadcaster. Be aware that the last thread has a
            ! different redSiz, so to calculate correctly we must use siz/nProc
            calcProc = (jj - 1)/(siz/nProc)
            IF (calcProc.GE.nProc) calcProc = nProc - 1
            
            ! Broadcast
            CALL MPI_BCAST(prevSol(:, jj), nShells, MPI_DOUBLE_PRECISION, &
                           calcProc, MPI_COMM_WORLD, ierror)
        END DO
    END IF
    
    ! Copy values to ovShell
    DO jj = 1, nShells
        ovShell(firstOv + jj - 1)%dens = prevSol(jj, :)
    END DO
    
    ! Update intShell values
    ii = 1
    entered = .FALSE.
    DO jj = 1, SIZE(ovShell)
        ! Copy convectives
        DO WHILE (intShell(ii)%isConvective)
            IF (.NOT.entered) entered = .TRUE.
            intShell(ii)%dens = ovShell(jj)%dens
            
            ii = ii + 1
            IF (ii.GE.totShell) EXIT
        END DO
        
        ! Correct index or copy non-convective
        IF (entered) THEN
            ii = ii - 1
            entered = .FALSE.
        ELSE
            intShell(ii)%dens = ovShell(jj)%dens
        END IF
        
        ! Advance index
        ii = ii + 1
    END DO
    
    ! Check if negatives exist and make them 0
    DO ii = 1, totShell - 1
        DO kk = 1, siz
            IF (ABS(intShell(ii)%dens(kk)).LT.yscale) THEN
                intShell(ii)%dens(kk) = 0.D0
            ELSE IF (intShell(ii)%dens(kk).LT.0.D0) THEN
                PRINT*, "# Warning, making zero shell and element ", ii, kk
                PRINT*, "# Previous value: ", intShell(ii)%dens(kk)
                intShell(ii)%dens(kk) = 0.D0
            ELSE IF (intShell(ii)%dens(kk).GT.1.D0) THEN
                PRINT*, "Big value after overshooting", firstOv
                PRINT*, ii, kk, intShell(ii)%dens(kk)
                STOP
            END IF
        END DO
    END DO
    
    ! Free memory
    DEALLOCATE(redSol)
    
END SUBROUTINE applyOvershooting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine solves the advective overshooting equations with a       !!!
!!! Richardson extrapolation method on an euler method.                      !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -ovMatrix, the overshooting matrix.                                      !!!
!!! -prevSol, initial and subsequent values.                                 !!!
!!! -nShells, size of ovMatrix.                                              !!!
!!! -dtot, total timestep.                                                   !!!
!!! -convecIndex, the array holding the indices for the convective zones as  !!!
!!!               well as the furthermost index they affect.                 !!!
!!! -yscale, threshold for zero abundances.                                  !!!
!!! -eps, relative error.                                                    !!!
!!! -siz, number of species.                                                 !!!
!!!                                                                          !!!
!!! On output prevSol holds the solution to the overshooting problem.        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solveAdvOv(ovMatrix, prevSol, nShells, dtot, convecIndex, yscale, &
                      eps, siz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::ovMatrix(:, :), prevSol(:, :), dtot, yscale, eps
    INTEGER::nShells, convecIndex(:, :), siz
    
    ! Local
    DOUBLE PRECISION::eta(5), numSol(nShells, siz, 5), ext(nShells, siz, 9)
    DOUBLE PRECISION::identity(nShells, nShells), systm(nShells, nShells)
    DOUBLE PRECISION::tMatx(nShells, nShells), bb(siz, nShells), error(siz)
    DOUBLE PRECISION::Hnew, stepErr, Htot, HH, xx(siz, nShells), errcoef
    DOUBLE PRECISION::k1(5), k2(5), nn(5), c1(5), normFactors(siz)
    DOUBLE PRECISION::transpOv(nShells, nShells)
    INTEGER::steps(5), ii, jj, sz, kk
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Define steps
    steps = (/1, 6, 10, 14, 22/)
    
    ! Define eta
    sz = SIZE(eta)
    DO ii = 1, sz
        eta(ii) = 1.D0/steps(ii)
    END DO
    
    ! Identity matrix
    identity = 0.D0
    DO ii = 1, nShells
        identity(ii, ii) = 1.D0
    END DO
    
    ! Calculate beforehand all coefficients
    k1 = 0; k2 = 0; nn = 0; c1 = 0
    c1(2:) = 1/(1 - eta(2:)/eta(1:sz - 1))
    
    k1(2:) = -eta(2:)*eta(1:sz - 1)
    
    k2(3:) = -k1(3:)*(eta(3:) + eta(2:sz - 1))
    k2(3:) = k2(3:) + k1(2:sz - 1)*(eta(2:sz - 1) + eta(1:sz - 2))
    k2(3:) = k2(3:)/(1 - k1(3:)/k1(2:sz - 1))
    
    nn(3:) = -(eta(2:sz - 1)**4 - eta(2:sz - 1)*eta(1:sz - 2)**3)*c1(2:sz - 1)
    nn(3:) = nn(3:) + (eta(3:)**4 - eta(3:)*eta(2:sz - 1)**3)*c1(3:)
    nn(3:) = nn(3:)/(1 - k1(3:)/k1(2:sz - 1))
    
    ! Error coefficient
    errcoef = nn(5) - nn(4)*k2(5)/k2(4)
    errcoef = ABS(errcoef/(nn(5) - nn(4)*(k2(5)/k2(4) + 1) + nn(3)*k2(4)/k2(3)))
    
    ! Get normalization factors and normalize
    DO ii = 1, siz
        normFactors(ii) = MAXVAL(prevSol(:, ii))
        IF (normFactors(ii).GT.0.D0) THEN
            prevSol(:, ii) = prevSol(:, ii)/normFactors(ii)
        END IF
    END DO
    
    ! Timesteps
    transpOv = TRANSPOSE(ovMatrix)
    Htot = dtot; HH = Htot
    DO WHILE (Htot.GT.0.D0)
        ! Solve system for each eta
        DO ii = 1, SIZE(steps)
            ! Define system
            systm = -HH*transpOv/steps(ii)
            DO jj = 1, nShells
                systm(jj, jj) = 1.D0 + systm(jj, jj)
            END DO
            
            ! Factorize it
            CALL sparseOvMatrixGauss(systm, tMatx, convecIndex, nShells)
            
            ! Solve system
            numSol(:, :, ii) = prevSol
            
            bb = TRANSPOSE(numSol(:, :, ii))
            DO jj = 1, steps(ii)
                ! Solve linear system systm*xx = bb
                CALL dimBackSub(systm, tMatx, bb, xx, convecIndex, nShells, siz)
                
                bb = xx
            END DO
            numSol(:, :, ii) = TRANSPOSE(xx)
        END DO
        
        ! Extrapolate
        
        ! First order extrapolations
        DO ii = 2, SIZE(numSol(1, 1, :))
            jj = ii - 1
            ext(:, :, jj) = eta(ii - 1)*numSol(:, :, ii)
            ext(:, :, jj) = ext(:, :, jj) - eta(ii)*numSol(:, :, ii - 1)
            ext(:, :, jj) = ext(:, :, jj)/(eta(ii - 1) - eta(ii))
        END DO
        
        ! Second order extrapolations
        DO ii = 3, SIZE(numSol(1, 1, :))
            jj = ii + SIZE(numSol(1, 1, :)) - 3
            kk = jj - SIZE(numSol(1, 1, :)) + 2
            ext(:, :, jj) = k1(ii - 1)*ext(:, :, kk) - k1(ii)*ext(:, :, kk - 1)
            ext(:, :, jj) = ext(:, :, jj)/(k1(ii - 1) - k1(ii))
        END DO
        
        ! Third order extrapolations
        DO ii = 4, SIZE(numSol(1, 1, :))
            jj = ii + 2*SIZE(numSol(1, 1, :)) - 6
            kk = jj - SIZE(numSol(1, 1, :)) + 3
            ext(:, :, jj) = k2(ii - 1)*ext(:, :, kk) - k2(ii)*ext(:, :, kk - 1)
            ext(:, :, jj) = ext(:, :, jj)/(k2(ii - 1) - k2(ii))
        END DO
        
        ! Calculate relative error
        DO jj = 1, siz
            error(jj) = 0.D0
            DO ii = 1, nShells
                stepErr = ext(ii, jj, 9) - ext(ii, jj, 8)
                
                IF (ABS(ext(ii, jj, 8)).GT.yscale) THEN
                    stepErr = stepErr/ext(ii, jj, 8)
                ELSE
                    stepErr = stepErr/yscale
                END IF
                
                error(jj) = error(jj) + stepErr**2
            END DO
        END DO
        error = SQRT(error/nShells)*errcoef
        
        ! Calculate new HH
        Hnew = HH*(0.25*eps/MAXVAL(error))**0.25
        
        ! If MAXVAL(error) > eps, start anew
        IF (MAXVAL(error).GT.eps) THEN
            HH = Hnew
            CYCLE
        END IF
        
        ! Else, advance solution
        prevSol(:, :) = ABS(ext(:, :, 9))
        
        ! Advance timestep
        Htot = Htot - HH
        HH = Hnew
        
        ! Correct timestep
        IF (HH.GT.Htot) HH = Htot
    END DO
    
    ! De-normalize
    DO ii = 1, siz
        IF (normFactors(ii).GT.0.D0) THEN
            prevSol(:, ii) = prevSol(:, ii)*normFactors(ii)
        END IF
    END DO
    
END SUBROUTINE solveAdvOv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This subroutine solves the diffusive overshooting equations with a       !!!
!!! Crank-Nicholson method.                                                  !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -ovMatrix, the overshooting matrix.                                      !!!
!!! -prevSol, initial and subsequent values.                                 !!!
!!! -nShells, size of ovMatrix.                                              !!!
!!! -dtot, total timestep.                                                   !!!
!!! -yscale, threshold for zero abundances.                                  !!!
!!! -eps, relative error.                                                    !!!
!!! -siz, number of species.                                                 !!!
!!!                                                                          !!!
!!! On output prevSol holds the solution to the overshooting problem.        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE solveDifOv(ovMatrix, prevSol, nShells, dtot, yscale, eps, siz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::ovMatrix(:, :), prevSol(:, :), dtot, yscale, eps
    INTEGER::nShells, siz
    
    ! Local
    DOUBLE PRECISION::Htot, HH, xx(nShells, siz), normFactors(siz), error
    DOUBLE PRECISION::tempErr, dtStep, coarseSol(nShells, siz)
    INTEGER::ii
    LOGICAL::firstTime
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Get normalization factors and normalize
    DO ii = 1, siz
        normFactors(ii) = MAXVAL(prevSol(:, ii))
        IF (normFactors(ii).GT.0.D0) THEN
            prevSol(:, ii) = prevSol(:, ii)/normFactors(ii)
        END IF
    END DO
    
    ! Define initial timestep
    firstTime = .TRUE.
    Htot = dtot; HH = dtot; dtStep = HH
    DO WHILE(Htot.GT.0)
        
        ! If too many steps, reduce the total timestep fast
        IF (HH/dtStep.GT.2.D2) THEN
            HH = HH*0.1
            dtStep = HH
            firstTime = .TRUE.
        END IF
        
        ! Solve system for HH
        xx = solveCrankThomas(ovMatrix, prevSol, nShells, HH, dtStep, siz)
        
        ! If first time, return now
        IF (firstTime) THEN
            firstTime = .FALSE.
            coarseSol = xx
            dtStep = dtStep*0.5
            CYCLE
        END IF
        
        ! Compare coarseSol and xx. If the difference is lower than eps,
        ! then this solution has converged
        
        ! Calculate error
        error = 0.D0
        DO ii = 1, siz
            tempErr = SUM(ABS((coarseSol(:, ii) - xx(:, ii))/(coarseSol(:, ii) &
                                                             + yscale)))/nShells
            IF (tempErr.GT.error) error = tempErr
        END DO
        
        ! Check if increase timestep, else reduce it
        IF (error.LE.eps) THEN
            prevSol = xx
            Htot = Htot - HH
            HH = MIN(HH*SQRT(eps/error), Htot)
            dtStep = HH
            firstTime = .TRUE.
        ELSE
            dtStep = dtStep*0.5
            coarseSol = xx
        END IF
    END DO
    
    ! De-normalize
    DO ii = 1, siz
        IF (normFactors(ii).GT.0.D0) THEN
            prevSol(:, ii) = prevSol(:, ii)*normFactors(ii)
        END IF
    END DO
    
END SUBROUTINE solveDifOv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This function solves the specific Crank-Nicholson implementation with a  !!!
!!! Thomas algorithm for a total time Htot and a timestep of HH.             !!!
!!!                                                                          !!!
!!! The input values are:                                                    !!!
!!! -ovMatrix, the overshooting matrix.                                      !!!
!!! -prevSol, initial and subsequent values.                                 !!!
!!! -nShells, size of ovMatrix.                                              !!!
!!! -dtot, total integration time.                                           !!!
!!! -dtStep, fractional timestep chosen.                                     !!!
!!! -siz, number of species.                                                 !!!
!!!                                                                          !!!
!!! On output prevSol holds the solution to the overshooting problem.        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION solveCrankThomas(ovMatrix, prevSol, nShells, dtot, dtStep, siz)
    IMPLICIT NONE
    
    ! Input
    DOUBLE PRECISION::ovMatrix(:, :), prevSol(:, :), dtot, dtStep
    INTEGER::nShells, siz
    
    ! Function
    DOUBLE PRECISION::solveCrankThomas(nShells, siz)
    
    ! Local
    DOUBLE PRECISION::xx(nShells, siz), coefs(nShells, 3), indep(nShells, siz)
    DOUBLE PRECISION::Htot, HH
    INTEGER::ii, jj
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF DECLARATIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    HH = dtStep
    Htot = dtot
    
    ! Calculate and solve first the coefficients matrix
    coefs(:, 1) = -HH*ovMatrix(:, 1)*0.5
    coefs(:, 2) = 1 + HH*ovMatrix(:, 2)*0.5
    coefs(:, 3) = -HH*ovMatrix(:, 3)*0.5
    
    xx = prevSol
    
    ! Now solve matrix with the Thomas algorithm
    coefs(1, 1) = coefs(1, 1)/coefs(1, 2)
    DO ii = 2, nShells - 1
        coefs(ii, 1) = coefs(ii, 1)/(coefs(ii, 2) - &
                                       coefs(ii, 3)*coefs(ii - 1, 1))
    END DO
    
    ! Finally solve diffusive step
    DO WHILE (Htot.GT.0.D0)
        ! Create independent array
        indep(1, :) = xx(1, :)*(1 - HH*ovMatrix(1, 2)*0.5) + &
                      xx(2, :)*HH*ovMatrix(1, 1)*0.5
        DO ii = 2, nShells - 1
            indep(ii, :) = xx(ii - 1, :)*HH*ovMatrix(ii, 3)*0.5 + &
                           xx(ii, :)*(1 - HH*ovMatrix(ii, 2)*0.5) + &
                           xx(ii + 1, :)*HH*ovMatrix(ii, 1)*0.5
        END DO
        indep(nShells, :) = xx(nShells, :)*(1 - HH*ovMatrix(nShells, 2)*0.5) + &
                            xx(nShells - 1, :)*HH*ovMatrix(nShells, 3)*0.5
        
        ! Now transform it
        DO ii = 1, siz
            indep(1, ii) = indep(1, ii)/coefs(1, 2)
            DO jj = 2, nShells
                indep(jj, ii) = indep(jj, ii) - &
                                coefs(jj, 3)*indep(jj - 1, ii)
                indep(jj, ii) = indep(jj, ii)/(coefs(jj, 2) - &
                                            coefs(jj, 3)*coefs(jj - 1, 1))
            END DO
        END DO
        
        ! And use it to solve
        DO ii = 1, siz
            xx(nShells, ii) = indep(nShells, ii)
            DO jj = nShells - 1, 1, -1
                xx(jj, ii) = indep(jj, ii) - xx(jj + 1, ii)*coefs(jj, 1)
            END DO
        END DO
        
        Htot = Htot - HH
        IF (HH.GT.Htot) THEN
            HH = Htot
            
            ! Calculate and solve the coefficients matrix with the new HH
            coefs(:, 1) = -HH*ovMatrix(:, 1)*0.5
            coefs(:, 2) = 1 + HH*ovMatrix(:, 2)*0.5
            coefs(:, 3) = -HH*ovMatrix(:, 3)*0.5
            
            ! Now solve matrix with the Thomas algorithm
            coefs(1, 1) = coefs(1, 1)/coefs(1, 2)
            DO ii = 2, nShells - 1
                coefs(ii, 1) = coefs(ii, 1)/(coefs(ii, 2) - &
                                               coefs(ii, 3)*coefs(ii - 1, 1))
            END DO
        END IF
    END DO
    
    solveCrankThomas = xx
    RETURN
    
END FUNCTION solveCrankThomas

END MODULE mixer

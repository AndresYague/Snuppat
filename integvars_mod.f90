!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                            !!
!! This module keeps variables for the integration method.                    !!
!!                                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE integvars_mod
    IMPLICIT NONE
    SAVE
    
    ! Substeps for the explicit method
    INTEGER, DIMENSION(18)::substeps = (/2, 3, 4, 6, 8, 12, 16, 24, 32, 48, &
                                         64, 96, 128, 192, 256, 384, 512, 768/)
    INTEGER, DIMENSION(18)::aj = 1
    DOUBLE PRECISION, DIMENSION(18, 18)::alph = 0
    INTEGER::nSubsteps = 16
    
END MODULE integvars_mod

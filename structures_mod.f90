!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                            !!
!! This module stores all the structures.                                     !!
!!                                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE structures_mod
    IMPLICIT NONE
    SAVE
    
    TYPE REACT
        DOUBLE PRECISION, ALLOCATABLE::avector(:)
        DOUBLE PRECISION, ALLOCATABLE::crossTable(:)
        INTEGER::targnum, prodnum, targindx(3), prodindx(4), totnum
        CHARACTER(10)::source
        LOGICAL::isEc
        TYPE (REACT), POINTER::next
    END TYPE REACT
    
    TYPE CROSSARR
        DOUBLE PRECISION::crossect
        LOGICAL::isEc, isRepeated
        INTEGER::targnum, prodnum, targindx(3), prodindx(4), totnum
        TYPE (CROSSARR), POINTER::next
    END TYPE CROSSARR
    
    TYPE SHELL
        DOUBLE PRECISION::temp, rho, mass, radiat, radius, hp, pressure
        DOUBLE PRECISION::velocity
        DOUBLE PRECISION, ALLOCATABLE::dens(:)
    END TYPE SHELL
    
    TYPE INTERSHELL
        TYPE (INTERSHELL), POINTER::convShell
        DOUBLE PRECISION::temp, rho, mass0, mass1, rad0, rad1, hp0, hp1, pres0
        DOUBLE PRECISION::pres1, vel0, vel1, radiat0, radiat1
        DOUBLE PRECISION, ALLOCATABLE::dens(:)
        INTEGER::ii0, ii1
        LOGICAL::isConvective = .FALSE.
    END TYPE INTERSHELL
    
    TYPE INDICES
        CHARACTER(5)::elem
        INTEGER::posit
    END TYPE INDICES
    
END MODULE structures_mod

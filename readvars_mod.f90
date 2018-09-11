!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                            !!
!! This module stores all the reading variables.                              !!
!!                                                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE readvars_mod
    IMPLICIT NONE
    SAVE
    
    ! Lists with croefficients (.dat and .lst are used for reactions_list and
    ! .lst are used for reactions).
    CHARACTER(3), DIMENSION(9), PARAMETER::filenames = (/"1a1", "1a2", "1a3", &
                                                         "2a1", "2a2", "2a3", &
                                                         "2a4", "3a1", "3a2"/)
    CHARACTER(10), DIMENSION(9), PARAMETER::lowTempFiles = (/"1a1lowTemp", &
                                            "1a2lowTemp", "1a3lowTemp", &
                                            "2a1lowTemp", "2a2lowTemp", &
                                            "2a3lowTemp", "2a4lowTemp", &
                                            "3a1lowTemp", "3a2lowTemp"/)
    
    ! Cross section tables temperature
    DOUBLE PRECISION, DIMENSION(60)::tempTable = 0.D0
    INTEGER::tabSiz = 60
    
    ! List with species and their labels.
    ! Physics and chemistry files
    CHARACTER(20), PARAMETER::species = "data/species.dat"
    CHARACTER(20), PARAMETER::atonSpecies = "data/species.dat13"
    CHARACTER(20), PARAMETER::physics = "input/physics.dat"
    CHARACTER(20), PARAMETER::chemistry = "input/chemistry.dat"
    CHARACTER(20), PARAMETER::output = "output/processed.dat"
    CHARACTER(30), PARAMETER::solar = "data/initial_compAGS2009.dat"
    CHARACTER(20), PARAMETER::contChem = "input/contChem.dat"
    
    ! Data reading utilities:
    INTEGER, PARAMETER::uni = 16, uni2 = 17, uni3 = 18
    INTEGER::error = 0
    
    ! List indication:
    INTEGER, DIMENSION(9), PARAMETER::targ = (/1, 1, 1, 2, 2, 2, 2, 3, 3/)
    INTEGER, DIMENSION(9), PARAMETER::prod = (/1, 2, 3, 1, 2, 3, 4, 1, 2/)

END MODULE readvars_mod

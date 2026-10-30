! GRASPG compile-time capacities and initial defaults; rebuild after editing.
! Shared radial/orbital capacities are defined only in suite_parameters.
module csfg_parameters
    use vast_kind_param, only: double, long
    implicit none
    ! NINCOR is a matrix-element count, not the interactive GB memory limit.
    integer(long), parameter :: CSFG_INCORE_MAX_ELEMENTS = 2147483648_long
    integer, parameter :: CSFG_MAX_BLOCKS = 50
    integer, parameter :: CSFG_INTERACT_MAX_BLOCKS = 20
    integer, parameter :: CSFG_MAX_SPLIT_LAYERS = 50
    integer, parameter :: CSFG_LABEL_LENGTH = 1070
    integer, parameter :: CSFG_MAX_SELECTED_LEVELS = 100
    integer, parameter :: CSFG_LEVEL_MAX_STATES = 20000
    integer, parameter :: CSFG_LEVEL_MAX_FILES = 1000
    ! JMax is tied to the JFraction lookup table in rlevels_csfg.
    integer, parameter :: CSFG_LEVEL_MAX_J = 22
    real(double), parameter :: CSFG_MATRIX_CUTOFF = 1.0d-10
    real(double), parameter :: CSFG_RESTART_CUTOFF = 1.0d-20
    logical, parameter :: CSFG_DEFAULT_LSYMD = .false.
    logical, parameter :: CSFG_DEFAULT_NDISCARDBR = .false.
    logical, parameter :: CSFG_DEFAULT_LDISCARDBR = .false.
    integer, parameter :: CSFG_DEFAULT_NDISBR = 15
    integer, parameter :: CSFG_DEFAULT_LDISBR = 14
    logical, parameter :: CSFG_DEFAULT_LCRITE = .false.
    real(double), parameter :: CSFG_DEFAULT_CRITE0 = 1.0d-12
end module csfg_parameters

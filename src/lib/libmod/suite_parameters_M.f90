! Shared compile-time limits. Keep upstream values as the default configuration.
! Edit here (or use the grid patcher), then rebuild all libraries and programs.
module suite_parameters
    implicit none
    private
    public :: KEYORB, NNNP, NNN1, NNNW, NNNWM1, NNNWM2

    integer, parameter :: KEYORB = 215
    integer, parameter :: NNNP = 590
    integer, parameter :: NNN1 = NNNP + 10
    integer, parameter :: NNNW = 127
    integer, parameter :: NNNWM1 = NNNW - 1
    integer, parameter :: NNNWM2 = NNNW - 2
end module suite_parameters

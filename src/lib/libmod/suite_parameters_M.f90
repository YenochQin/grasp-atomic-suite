! Shared compile-time limits. Keep upstream values as the default configuration.
! Edit here (or use the grid patcher), then rebuild all libraries and programs.
module suite_parameters
    implicit none
    private
    public :: keyorb, nnnp, nnn1, nnnw, nnnwm1, nnnwm2

    integer, parameter :: keyorb = 215
    integer, parameter :: nnnp = 590
    integer, parameter :: nnn1 = nnnp + 10
    integer, parameter :: nnnw = 127
    integer, parameter :: nnnwm1 = nnnw - 1
    integer, parameter :: nnnwm2 = nnnw - 2
end module suite_parameters

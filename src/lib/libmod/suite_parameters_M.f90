! The single source of shared compile-time limits and numerical defaults.
! Edit this file, then rebuild ALL libraries and serial/MPI applications.
module suite_parameters
    use vast_kind_param, only: DOUBLE
    implicit none
    private
    public :: KEYORB, NNNP, NNN1, NNNW, NNNWM1, NNNWM2
    public :: FINITE_RNT_SCALE, FINITE_H, FINITE_N
    public :: POINT_RNT_SCALE, POINT_H, POINT_N, DEFAULT_HP, DEFAULT_ACCY, NODE_THRESHOLD

    ! Array capacities and packing limits. Derived values must not be edited.
    integer, parameter :: KEYORB = 215
    integer, parameter :: NNNP = 2990
    integer, parameter :: NNN1 = NNNP + 10
    integer, parameter :: NNNW = 127
    integer, parameter :: NNNWM1 = NNNW - 1
    integer, parameter :: NNNWM2 = NNNW - 2

    ! RNT is the coefficient below divided by the nuclear charge Z.
    real(DOUBLE), parameter :: FINITE_RNT_SCALE = 2.0D-6
    real(DOUBLE), parameter :: FINITE_H = 5.0D-2
    integer, parameter :: FINITE_N = NNNP
    real(DOUBLE), parameter :: POINT_RNT_SCALE = EXP(-65.0D0/16.0D0)
    real(DOUBLE), parameter :: POINT_H = 0.0625D0
    integer, parameter :: POINT_N = MIN(220,NNNP)
    real(DOUBLE), parameter :: DEFAULT_HP = 0.0D0

    ! Zero selects adaptive H**6. A positive value sets a fixed default ACCY.
    ! Interactive ACCY and existing GRASP_COUNT_ACCY overrides retain priority.
    real(DOUBLE), parameter :: DEFAULT_ACCY = 0.0D0
    real(DOUBLE), parameter :: NODE_THRESHOLD = 0.05D0
end module suite_parameters

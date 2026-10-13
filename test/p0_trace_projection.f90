program p0_trace_projection
    use vast_kind_param, only: double
    use parameter_def, only: nnnp
    use orb_C, only: np, nh, nak, e, gama
    use orba_C, only: iorder
    use corre_C, only: lcorre
    use fixd_C, only: lfix
    use invt_C, only: noinvt
    use wave_C, only: pf, qf, pz, mf
    use scf_C, only: method, scnsty
    use grid_C, only: r
    use def_C, only: accy
    use p0_trace
    use orthy_I
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    character(len=1024) :: directory
    character(len=16) :: argument
    integer :: rank
    real(double) :: expected
    call get_command_argument(1, directory)
    call get_command_argument(2, argument)
    read(argument,*) rank
    allocate(pf(nnnp,2), qf(nnnp,2))
    pf = 0.d0
    qf = 0.d0
    pf(2,1) = 1.d0
    pf(2,2) = 0.5d0
    pf(3,2) = sqrt(0.75d0)
    pz(:2) = [1.d0, 2.d0]
    mf(:2) = 3
    np(:2) = [1,2]
    nh(:2) = 's '
    nak(:2) = -1
    e(:2) = [2.d0, 1.d0]
    gama(:2) = 1.d0
    iorder(:2) = [1,2]
    lfix(:2) = [.true., .false.]
    lcorre(:2) = [.false., .true.]
    noinvt(:2) = .true.
    method(:2) = [1,3]
    scnsty(:2) = [0.d0, 0.1d0]
    r = 0.d0
    r(2:7) = [1.d0,2.d0,3.d0,4.d0,5.d0,6.d0]
    accy = 1.d-8
    call p0_trace_open(trim(directory), rank)
    call p0_trace_cycle(2)
    call p0_trace_update(2)
    call p0_trace_iteration(3)
    call p0_trace_trial('trial_probe', 2, 2.d0, pf(:,2), qf(:,2), 3, alpha=0.25d0)
    if (pz(2) /= 2.d0 .or. pf(2,2) /= 0.5d0) error stop 'trace changed input'
    ! Production ORTHY, with the exact two-dimensional integral seam below.
    call orthy(2, 2, .false.)
    expected = 1.5d0 / sqrt(0.75d0)
    if (abs(pz(2)-expected) > 1.d-14) error stop 'incorrect coefficient after projection'
    if (pf(2,2) /= 0.d0 .or. abs(pf(3,2)-1.d0) > 1.d-14) error stop 'incorrect projection'
    if (pz(1) /= 1.d0 .or. pf(2,1) /= 1.d0) error stop 'fixed orbital changed'
    if (any(qf /= 0.d0)) error stop 'small component changed'
    pz(2) = ieee_value(0.d0, ieee_quiet_nan)
    call p0_trace_stored('nonfinite_probe', 2)
    if (.not.ieee_is_nan(pz(2))) error stop 'trace replaced nonfinite value'
    call p0_trace_close()
    print *, 'Projection and trace-neutrality checks passed.'
end program p0_trace_projection

real(kind(0.d0)) function rint(i, j, k)
    use wave_C, only: pf, qf
    implicit none
    integer, intent(in) :: i, j, k
    if (k /= 0) error stop 'unexpected integral in projection fixture'
    rint = sum(pf(2:3,i)*pf(2:3,j) + qf(2:3,i)*qf(2:3,j))
end function rint

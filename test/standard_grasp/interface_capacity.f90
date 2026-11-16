! The compiler must reject procedure association if an imported interface still
! contains a private radial/orbital extent instead of the suite configuration.
program standard_interface_capacity
    use parameter_def, only: NNNP, NNN1, NNNW
    use, intrinsic :: iso_fortran_env, only: real64
    use tail_I, only: tail
    use frmtfp_I, only: frmtfp
    use frmrwf_I, only: frmrwf
    implicit none (type, external)
    procedure(tail), pointer :: radial
    procedure(frmtfp), pointer :: estimate
    procedure(frmrwf), pointer :: imported
    real(real64) :: p(NNNP), q(NNNP)
    integer :: indices(NNNW), mtp

    if (NNNP /= 790 .or. NNN1 /= 800 .or. NNNW /= 149) error stop 'wrong test configuration'
    radial => probe_tail
    estimate => probe_estimate
    imported => probe_import
    p = 0
    indices = 0
    call radial(1, p, q, 1, mtp)
    call estimate(indices, NNNW)
    call imported(indices, NNNW, 'fixture')
    if (mtp /= NNNP .or. q(NNNP) /= 1) error stop 'radial interface mismatch'
contains
    subroutine probe_tail(iorb, p, q, jp, mtp)
        integer, intent(in) :: iorb, jp
        integer, intent(out) :: mtp
        real(real64), intent(inout) :: p(NNNP)
        real(real64), intent(out) :: q(NNNP)
        q = 0
        q(NNNP) = 1
        mtp = NNNP
    end subroutine probe_tail

    subroutine probe_estimate(index, nsubs)
        integer, intent(in) :: index(NNNW), nsubs
        if (nsubs /= size(index)) error stop 'estimate interface mismatch'
    end subroutine probe_estimate

    subroutine probe_import(index, nsubs, filnam)
        integer, intent(in) :: index(NNNW), nsubs
        character(*), intent(in) :: filnam
        if (nsubs /= size(index) .or. filnam /= 'fixture') error stop 'import interface mismatch'
    end subroutine probe_import
end program standard_interface_capacity

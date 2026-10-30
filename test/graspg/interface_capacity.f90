! Procedure pointer assignment checks explicit interface shapes at compile time.
! A private copy of NNNP/NNNW in either interface breaks the alternate-capacity test.
program csfg_interface_capacity
    use parameter_def, only: NNNP, NNNW
    use vast_kind_param, only: double
    use tail_I, only: tail
    use frmtfp_I, only: frmtfp
    implicit none
    procedure(tail), pointer :: radial
    procedure(frmtfp), pointer :: orbital
    real(double) :: p(NNNP), q(NNNP)
    integer :: indices(NNNW), mtp
    radial => probe_tail
    orbital => probe_frmtfp
    p = 0
    indices = 0
    call radial(1, p, q, 1, mtp)
    call orbital(indices, NNNW)
    if (mtp /= NNNP .or. q(NNNP) /= 1) error stop 'radial interface mismatch'
contains
    subroutine probe_tail(iorb, p, q, jp, mtp)
        integer, intent(in) :: iorb, jp
        integer, intent(out) :: mtp
        real(double), intent(inout) :: p(NNNP)
        real(double), intent(out) :: q(NNNP)
        q = 0
        q(NNNP) = 1
        mtp = NNNP
    end subroutine
    subroutine probe_frmtfp(index, nsubs)
        integer, intent(in) :: index(NNNW), nsubs
        if (nsubs /= size(index)) error stop 'orbital interface mismatch'
    end subroutine
end program csfg_interface_capacity

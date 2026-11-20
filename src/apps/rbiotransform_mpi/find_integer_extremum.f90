!
!     ------------------------------------------------------------------
!     I F N M N X
!     ------------------------------------------------------------------
!
integer function find_integer_extremum (ivec, nel, imxmn)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !
    ! Smallest or largest value of integer array
    !
    ! IMXMN = 1 => Largest value
    ! IMXMN = 2 => Smallest
    !
    !************************************************************************
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: nel
    integer , intent(in) :: imxmn
    integer , intent(in) :: ivec(nel)

    ! Local declarations
    integer :: ival, iel, ntest
    !-----------------------------------------------
    !
    ival = ivec(1)
    if (imxmn == 1) then
        ival = max0(maxval(ivec(2:nel)),ival)
    else if (imxmn == 2) then
        ival = min0(minval(ivec(2:nel)),ival)
    else
        write (output_unit, *) " Stop in IFNMNX "
        write (output_unit, *) " Improper calue of IMXMN ", imxmn
        stop "IFNMNX"
    end if
    !
    find_integer_extremum = ival
    !
    ntest = 0
    if (ntest /= 0) write (output_unit, *) " Value returned from IFNMNX", find_integer_extremum
    !
    return
end function find_integer_extremum

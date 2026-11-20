!
!     ------------------------------------------------------------------
!     I E L S U M
!     ------------------------------------------------------------------
!
integer function sum_integer_vector (ivec, nelmnt)
    !
    ! Sum elements of integer array
    !
    !************************************************************************
    !-----------------------------------------------
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: nelmnt
    integer, intent(in) :: ivec(*)

    ! Local declarations
    integer :: isum, iel
    !-----------------------------------------------
    !
    isum = 0
    isum = sum(ivec(:nelmnt))
    !
    sum_integer_vector = isum
    !
    return
end function sum_integer_vector

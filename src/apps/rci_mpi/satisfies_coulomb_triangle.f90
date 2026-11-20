!***********************************************************************
!                                                                      *
logical function satisfies_coulomb_triangle (la, k, lb)
    !                                                                      *
    !   Written by Per Jonsson                                             *
    !                                                                      *
    !***********************************************************************
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: la
    integer, intent(in) :: k
    integer, intent(in) :: lb
    !
    !   Perform the triangularity check
    !
    if (mod(k + la + lb,2) /= 0) then
        satisfies_coulomb_triangle = .false.
    else
        if (abs(la - lb) > k) then
            satisfies_coulomb_triangle = .false.
        else if (la + lb < k) then
            satisfies_coulomb_triangle = .false.
        else
            satisfies_coulomb_triangle = .true.
        end if
    end if

    return
end function satisfies_coulomb_triangle

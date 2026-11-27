module transition_off_diagonal_flags
    implicit none (type, external)

    ! Module state
    integer :: m1_different_j_only_flag = 0, e2_different_j_only_flag = 0
contains
    logical function different_j_only(multipole)
        integer, intent(in) :: multipole
        different_j_only = (multipole == 1 .and. m1_different_j_only_flag == 1) .or. &
            (multipole == 2 .and. e2_different_j_only_flag == 1)
    end function different_j_only
end module transition_off_diagonal_flags

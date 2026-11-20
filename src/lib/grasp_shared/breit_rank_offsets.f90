! Stores list offsets by multipole rank for the two Breit integral groups.
module breit_rank_offsets
    implicit none (type, external)

    ! Module constants
    integer, parameter :: maximum_breit_multipole_rank = 20

    ! Module state
    integer :: first_breit_integral_group_offsets(0:maximum_breit_multipole_rank), &
        second_breit_integral_group_offsets(0:maximum_breit_multipole_rank)
end module breit_rank_offsets

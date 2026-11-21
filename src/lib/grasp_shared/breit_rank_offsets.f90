! Stores list offsets by multipole rank for the two Breit integral groups.
module breit_rank_offsets
    use grasp_parameters, only: maximum_orbital_two_j
    implicit none (type, external)

    ! Module constants
    integer, parameter :: maximum_breit_multipole_rank = maximum_orbital_two_j

    ! Module state
    integer :: first_breit_integral_group_offsets(0:maximum_breit_multipole_rank+1), &
        second_breit_integral_group_offsets(0:maximum_breit_multipole_rank+1)
end module breit_rank_offsets

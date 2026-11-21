module coulomb_rank_offsets
    use grasp_parameters, only: maximum_orbital_two_j
    implicit none (type, external)

    ! Module constants
    integer, parameter :: maximum_coulomb_multipole_rank = maximum_orbital_two_j

    ! Module state
    integer :: coulomb_integral_offsets_by_rank(0:maximum_coulomb_multipole_rank+1)
end module coulomb_rank_offsets

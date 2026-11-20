module coulomb_rank_offsets
    implicit none (type, external)

    ! Module constants
    integer, parameter :: maximum_coulomb_multipole_rank = 20

    ! Module state
    integer :: coulomb_integral_offsets_by_rank(0:maximum_coulomb_multipole_rank)
end module coulomb_rank_offsets

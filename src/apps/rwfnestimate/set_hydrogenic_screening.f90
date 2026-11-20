

!***********************************************************************
subroutine set_hydrogenic_screening(core_subshell_count)
    ! Module dependencies
    !
    ! Purpose:
    !   Compute hydrogenic screen parameters
    !
    ! Input:
    !   ncore
    !   in the common - nw, nkj()
    !
    ! Output:
    !   in the common - sigma()
    !
    !***********************************************************************
    use orbital_csf_basis, only: orbital_count, orbital_two_j
    use hydrogenic_screening, only: sigma
    implicit none (type, external)

    ! Arguments
    integer , intent(in) :: core_subshell_count

    ! Local declarations
    integer :: nelectron, i
    !-----------------------------------------------

    !...Core orbitals
    nelectron = 0
    do i = 1, core_subshell_count
        sigma(i) = nelectron + (orbital_two_j(i)+1)/2
        nelectron = nelectron + orbital_two_j(i) + 1
    end do

    !...Peel orbitals
    sigma(core_subshell_count+1:orbital_count) = nelectron

    return
end subroutine set_hydrogenic_screening

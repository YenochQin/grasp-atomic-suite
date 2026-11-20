module set_hydrogenic_screening_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: set_hydrogenic_screening
contains

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
        use hydrogenic_screening, only: orbital_screening_constants
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: core_subshell_count

        ! Local declarations
        integer :: nelectron, i
        !-----------------------------------------------

        !...Core orbitals
        nelectron = 0
        do i = 1, core_subshell_count
            orbital_screening_constants(i) = nelectron + (orbital_two_j(i)+1)/2
            nelectron = nelectron + orbital_two_j(i) + 1
        end do

        !...Peel orbitals
        orbital_screening_constants(core_subshell_count+1:orbital_count) = nelectron

        return
    end subroutine set_hydrogenic_screening
end module set_hydrogenic_screening_interface

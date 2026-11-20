module map_merged_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: map_merged_orbitals
contains

    !***********************************************************************
    !                                                                      *
    subroutine map_merged_orbitals
        ! Module dependencies
        !                                                                      *
        !   The position of an orbital in the merged list is connected to      *
        !   the positions in the initial and final state lists                 *
        !                                                                      *
        !   Written by Per Jonsson                                             *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use biorthogonal_orbital_sets, only: final_orbital_count, initial_orbital_count, final_orbital_principal_n, &
            final_orbital_kappa, initial_orbital_principal_n, initial_orbital_kappa, merged_to_final_orbital_index, &
            merged_to_initial_orbital_index
        use orbital_csf_basis, only: orbital_count, orbital_kappa, orbital_principal_n
        implicit none (type, external)

        ! Local declarations
        integer :: orbital_index, j
        !-----------------------------------------------
        !
        !ww      INTEGER PNTRIQ
        !     POINTER (PNTRIQ,RIQDUMMY)
        !
        !
        !     Initialize
        !
        merged_to_initial_orbital_index(:orbital_count) = 0.0_real64
        merged_to_final_orbital_index(:orbital_count) = 0.0_real64
        !
        !   Loop over the orbitals in the merged list
        !
        do orbital_index = 1, orbital_count
            do j = 1, initial_orbital_count
                if (orbital_principal_n(orbital_index)/=initial_orbital_principal_n(j) .or. &
                    orbital_kappa(orbital_index)/=initial_orbital_kappa(j)) cycle
                merged_to_initial_orbital_index(orbital_index) = j
            end do
            do j = 1, final_orbital_count
                if (orbital_principal_n(orbital_index)/=final_orbital_principal_n(j) .or. &
                    orbital_kappa(orbital_index)/=final_orbital_kappa(j)) cycle
                merged_to_final_orbital_index(orbital_index) = j
            end do
        end do

        return
    end subroutine map_merged_orbitals
end module map_merged_orbitals_interface

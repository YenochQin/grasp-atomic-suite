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
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use orbital_csf_basis, only: max_orbitals, orbital_count, orbital_kappa, orbital_principal_n
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
    nnii(:orbital_count) = 0.0_real64
    nnff(:orbital_count) = 0.0_real64
    !
    !   Loop over the orbitals in the merged list
    !
    do orbital_index = 1, orbital_count
        do j = 1, nwii
            if (orbital_principal_n(orbital_index)/=npii(j) .or. orbital_kappa(orbital_index)/=nakii(j)) cycle
            nnii(orbital_index) = j
        end do
        do j = 1, nwff
            if (orbital_principal_n(orbital_index)/=npff(j) .or. orbital_kappa(orbital_index)/=nakff(j)) cycle
            nnff(orbital_index) = j
        end do
    end do

    return
end subroutine map_merged_orbitals

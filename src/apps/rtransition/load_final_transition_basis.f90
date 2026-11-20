!***********************************************************************
!                                                                      *
subroutine load_final_transition_basis(final_core_subshell_count, state_name)
    ! Module dependencies
    !                                                                      *
    !   Open, check, load data from and close the  .csl  file. This file   *
    !   is always attached to stream 21.                                   *
    !                                                                      *
    !   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use atomic_calculation_state, only: electron_count, nelecff
    use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
        orbital_symmetry_labels, orbital_two_j
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii

    use csf_full_basis_reader, only: load_csf_basis
    use transition_basis_blocks, only: nblock, ncfblk, nblocki, nblockf, ncfi, ncff
    implicit none (type, external)

    ! Arguments
    integer , intent(out) :: final_core_subshell_count
    character , intent(in) :: state_name*24

    ! Local declarations
    integer :: j, io_status, ncorer, i
    character :: record*15
    !-----------------------------------------------
    !
    !
    !     Common relevant for the radial wave functions of the final state.
    !
    !
    !   The  .csl  file is FORMATTED; it must exist
    !
    call load_csf_basis(trim(state_name),ncorer,io_status)
    if (io_status /= 0) error stop "Transition descriptor load failed"
    if (nblock > size(ncff)) error stop "Transition block capacity exceeded"
    nblockf=nblock
    ncff=0
    ncff(:nblock)=ncfblk(:nblock)
    !
    !   Check if the core should be redefined
    !
    final_core_subshell_count = ncorer
    !      DO 3 I = NCORER+1,NW
    !         IFULLI = NKJ(I)+1
    !         DO 2 J = 1,NCF
    !    2    CONTINUE
    !         CALL CONVRT (NP(I),RECORD,LENTH)
    !         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
    !         PRINT *, ' in all CSFs; including this'
    !         PRINT *, ' subshell in the core;'
    !         NCORE = NCORE+1
    !    3 CONTINUE
    !
    nelecff = electron_count
    nwff = orbital_count
    ncfff = csf_count
    nhff(:orbital_count) = orbital_symmetry_labels(:orbital_count)
    npff(:orbital_count) = orbital_principal_n(:orbital_count)
    nakff(:orbital_count) = orbital_kappa(:orbital_count)
    nklff(:orbital_count) = orbital_l(:orbital_count)
    nkjff(:orbital_count) = orbital_two_j(:orbital_count)
    return
end subroutine load_final_transition_basis

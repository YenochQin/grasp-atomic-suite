!***********************************************************************
!                                                                      *
subroutine load_initial_transition_basis (initial_core_subshell_count,state_name)
    ! Module dependencies
    !                                                                      *
    !   Open, check, load data from and close the  .csl  file. This file   *
    !   is always attached to stream 21.                                   *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
    !                                                                      *
    !                                                                      *
    !***********************************************************************
    !...Translated by Charlotte Froese Fischer
    !                       Gediminas Gaigalas  10/05/17
    use grasp_parameters,   only: max_orbitals
    use array_memory_management, only: alloc, dalloc
    use atomic_calculation_state, only: electron_count, nelecii, nelecr
    use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
        nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, iqar, max_orbitals, nakr, ncfr, nhr, npr, &
        nwr, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, orbital_symmetry_labels, orbital_two_j
    use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
    use csf_quantum_queries, only: iq
    use csf_full_basis_reader, only: load_csf_basis
    use transition_basis_blocks, only: nblock, ncfblk, nblocki, nblockf, ncfi, ncff
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    integer :: initial_core_subshell_count
    character(len=24) :: state_name

    ! Local declarations
    integer :: i, j, io_status, core_subshell_count
    character(len=15) :: record
    !
    !   Check the first record of the file; if not as expected, try again
    !
    call load_csf_basis(trim(state_name),core_subshell_count,io_status)
    if (io_status /= 0) error stop "Transition descriptor load failed"
    if (nblock > size(ncfi)) error stop "Transition block capacity exceeded"
    nblocki=nblock
    ncfi=0
    ncfi(:nblock)=ncfblk(:nblock)
    !
    !   Check if the core should be redefined
    !
    initial_core_subshell_count = core_subshell_count
    !
    !   Copy data to reference arrays
    !
    nelecr = electron_count
    nelecii = electron_count
    ncfr = csf_count
    nwr = orbital_count
    !
    call alloc (iqar, max_orbitals, csf_count, "IQAR", "LDCSL1" )
    call alloc (jqsar,max_orbitals,3,csf_count, "JQSAR", "LDCSL1")
    call alloc (jcupar,max_orbitals, csf_count, "JCUPAR", "LDCSL1")
    !
    do i = 1,csf_count
        do j = 1,max_orbitals
            iqar(j,i) = csf_subshell_occupations(j,i)
            jqsar(j,1,i) = csf_subshell_quantum_codes(j,1,i)
            jqsar(j,2,i) = csf_subshell_quantum_codes(j,2,i)
            jqsar(j,3,i) = csf_subshell_quantum_codes(j,3,i)
            jcupar(j,i) = csf_coupling_quantum_codes(j,i)
        end do
    end do
    !
    nwii=orbital_count
    ncfii=csf_count
    do i = 1,orbital_count
        nhii(i) = orbital_symmetry_labels(i)
        npii(i) = orbital_principal_n(i)
        nakii(i) = orbital_kappa(i)
        nklii(i) = orbital_l(i)
        nkjii(i) = orbital_two_j(i)
        npr(i) = orbital_principal_n(i)
        nakr(i) = orbital_kappa(i)
        nhr(i) = orbital_symmetry_labels(i)
    end do
    !
    !   Deallocate storage so that LODCSL can restart correctly
    !
    call dalloc (csf_subshell_occupations, "IQA", "LDCSL1")
    call dalloc (csf_subshell_quantum_codes, "JQSA", "LDCSL1")
    call dalloc (csf_coupling_quantum_codes, "JCUPA", "LDCSL1")
    !


    return
end subroutine load_initial_transition_basis

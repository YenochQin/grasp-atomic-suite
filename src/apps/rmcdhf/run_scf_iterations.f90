!***********************************************************************
!                                                                      *
subroutine run_scf_iterations(extended_optimal_level, output_wavefunction_path)
    ! Module dependencies
    use calculation_file_units, only: mixing_output_unit
    !                                                                      *
    !   This  subroutine  performs  the SCF iterations. The procedure is   *
    !   essentially algorithm 5.1 of C Froese Fischer, Comput Phys Rep 3   *
    !   (1986) 290.                                                        *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
    !               [RSCF92]: improv, matrix, MAXARR, newco,               *
    !                         ORBOUT, ORTHSC, setlag.                      *
    !                                                                      *
    !   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
    !   Block version by Xinghong He          Last revision: 05 Aug 1998   *
    !   Midified by G. Gaigalas                              05 Feb 2017   *
    !      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
    !                                  JCUPA(NNNW*NCF)                     *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use array_memory_management, only: alloc, dalloc
    use transition_block_indices, only: idxblk
    use calculation_input_defaults, only: inpci, ndef, ndump
    use atomic_calculation_state, only: accy, iccmin, ncmin, nscf, nsic, wt
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    use fixed_orbital_flags, only: fixed_orbital_count, is_orbital_fixed, max_orbitals
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use orbital_lagrange_multipliers, only: ecv, iecc, nec
    use mcp_file_metadata, only: max_mcp_file_multipole_rank
    use serial_mpi_context, only: mpi_rank
    use symmetry_block_offsets, only: total_csf_count, mixing_coefficient_storage_size, &
        csf_offsets_by_block, eigenvalue_offsets_by_block, mixing_offsets_by_block
    use block_mean_energies, only: eavblk
    use orbital_csf_basis, only: csf_count, csf_subshell_occupations, max_orbitals, orbital_count
    use scf_orbital_order, only: iorder, max_orbitals
    use scf_iteration_state, only: max_orbitals, ucf, method, scnsty, epsmin, epsmax, emin, emax, zinf, ndcof, nxcof, &
        nycof, nddim, nxdim, nydim, da, xa, ya, nda, nxa, nya
    use scf_orthogonality_control, only: orthst

    ! Procedure interfaces
    use solve_hamiltonian_blocks_interface, only: solve_hamiltonian_blocks
    use update_scf_mixing_interface, only: update_scf_mixing
    use calculate_lagrange_multipliers_interface, only: calculate_lagrange_multipliers
    use improve_orbital_interface, only: improve_orbital
    use find_least_consistent_orbital_interface, only: find_least_consistent_orbital
    use print_radial_wavefunction_interface, only: print_radial_wavefunction

    ! Module dependencies
    use radial_orbital_orthogonalization, only: orthsc

    ! Procedure interfaces
    use write_scf_radial_wavefunctions_interface, only: write_scf_radial_wavefunctions
    use finish_summary_interface, only: finish_summary
    implicit none (type, external)

    ! Arguments
    !-----------------------------------------------
    !   G l o b a l   P a r a m e t e r s
    !-----------------------------------------------
    logical, intent(in) :: extended_optimal_level
    character  :: output_wavefunction_path*(*)

    ! Local declarations
    integer :: orbital_index, i, nit, jseq, kount, k
    real(real64) :: wtaev, wtaev0, dampmx
    logical :: convg, lsort, dvdfirst
    !-----------------------------------------------


    total_csf_count = csf_count
    !IF (myid .EQ. 0) PRINT *, '===SCF==='

    !=======================================================================
    !   Determine Orthonomalization order --- lsort
    !=======================================================================

    if (ndef == 0) then
        lsort = .false.
    else
        do
            write (stderr_unit, *) "Orthonomalization order? "
            write (stderr_unit, *) "     1 -- Update order"
            write (stderr_unit, *) "     2 -- Self consistency connected"
            read (stdin_unit, *) orbital_index
            if (orbital_index == 1) then
                lsort = .false.
            else if (orbital_index == 2) then
                lsort = .true.
            else
                write (stderr_unit, *) "Input is wrong, redo..."
                cycle
            end if
            exit
        end do
    end if

    !=======================================================================
    !   Deallocate storage that will no longer be used
    !=======================================================================

    !GG      CALL DALLOC (JQSA, 'JQSA', 'SCF')

    !=======================================================================
    !   Allocate and fill in auxiliary arrays
    !=======================================================================

    call alloc (csf_offsets_by_block, symmetry_block_count, "NCFPAST", "SCF")
    call alloc (eigenvalue_offsets_by_block, symmetry_block_count, "NCMINPAST", "SCF")
    call alloc (mixing_offsets_by_block, symmetry_block_count, "NEVECPAST", "SCF")
    call alloc (eavblk, symmetry_block_count, "EAVBLK", "SCF")

    csf_offsets_by_block(1) = 0
    eigenvalue_offsets_by_block(1) = 0
    mixing_offsets_by_block(1) = 0
    do i = 2, symmetry_block_count
        csf_offsets_by_block(i) = csf_offsets_by_block(i-1) + csf_counts_by_block(i - 1)
        eigenvalue_offsets_by_block(i) = eigenvalue_offsets_by_block(i-1) + eigenvector_counts_by_block(i - 1)
        mixing_offsets_by_block(i) = mixing_offsets_by_block(i-1) + eigenvector_counts_by_block(i - 1)*csf_counts_by_block(i - 1)
    end do

    !*** Size of the eigenvector array for all blocks
    mixing_coefficient_storage_size = mixing_offsets_by_block(symmetry_block_count) + &
        & eigenvector_counts_by_block(symmetry_block_count)*csf_counts_by_block(symmetry_block_count)

    if (extended_optimal_level) then
        call alloc (shifted_eigenvalues, ncmin, "EVAL", "SCF")
        call alloc (csf_mixing_coefficients, mixing_coefficient_storage_size, "EVEC", "SCF")
        !GG         CALL ALLOC (IATJPO, NCMIN, 'IATJPO', 'SCF')
        !GG         CALL ALLOC (IASPAR, NCMIN, 'IASPAR', 'SCF')
    end if

    !=======================================================================
    !
    !=======================================================================
    nddim = 0
    nxdim = 0
    nydim = 0

    !     This call should only be made AFTER the call to newco
    !     CALL setlag (EOL)

    !   For (E)OL calculations, determine the level energies and
    !   mixing coefficients

    !CFF   .. set the logical variable dvdfirst
    dvdfirst = .true.
    if (extended_optimal_level) then
        call solve_hamiltonian_blocks (dvdfirst)
        call update_scf_mixing (wtaev)
    end if
    wtaev0 = 0.0
    dvdfirst = .false.
    do nit = 1, nscf
        if (mpi_rank == 0) write (*, '(/," Iteration number ",1i3,/," --------------------")') nit

        !   For all pairs constrained through a Lagrange multiplier, compute
        !   the Lagrange multiplier

        call calculate_lagrange_multipliers (extended_optimal_level)

        !   Improve all orbitals in turn

        dampmx = 0.0
        if (mpi_rank == 0) write (*, '(41x,"Self-            Damping"/,"Subshell    Energy    Method   P0    ",&
            &"consistency  Norm-1  factor  JP"," MTP INV NNP"/)')
        do orbital_index = 1, orbital_count
            jseq = iorder(orbital_index)
            if (is_orbital_fixed(jseq)) cycle
            call improve_orbital (extended_optimal_level, jseq, lsort, dampmx)
        end do
        !
        !   For KOUNT = 1 to NSIC: find the least self-consistent orbital;
        !   improve it
        !
        !         write(istde,*) 'nsic=',nsic
        do kount = 1, nsic
            call find_least_consistent_orbital (k)
            if (k == 0) then
                convg = .true.
                exit
            else
                if (scnsty(k) <= accy) then
                    convg = .true.
                    exit
                end if
            end if
            call improve_orbital (extended_optimal_level, k, lsort, dampmx)
        end do

        if (kount > nsic) then
            call find_least_consistent_orbital (k)

            if (k == 0) then
                convg = .true.
            else
                if (scnsty(k) <= accy) then
                    convg = .true.
                else
                    convg = .false.
                end if
            end if

        end if
        if (ldbpr(24) .and. mpi_rank==0) call print_radial_wavefunction (0)

        !   Perform Gram-Schmidt process
        !   For OL calculation, orthst is true and orbitals are orthonormalized
        !   in subroutine improv. For AL calculation, orthst is false.
        if (.not.orthst) call orthsc

        !   Write the subshell radial wavefunctions to the .rwf file

        if (mpi_rank == 0) call write_scf_radial_wavefunctions (output_wavefunction_path)

        if (extended_optimal_level) then
            call solve_hamiltonian_blocks(dvdfirst)
            call update_scf_mixing (wtaev)
        end if
        !        Make this a relative convergence test
        !        IF(ABS(WTAEV-WTAEV0).LT.1.0D-9.and.
        !    &   DAMPMX.LT.1.0D-4) CONVG=.true.
        !        PRINT *, 'WTAEV, WTAEV0', WTAEV, WTAEV0,
        !    &             ABS((WTAEV-WTAEV0)/WTAEV)
        !cjb unified convergence criteria in RMCDHF and RMCDHF_MPI
        !cjb     IF(DABS(WTAEV-WTAEV0).LT.1.0D-8.and.              &
        !cjb                   DAMPMX.LT.1.0D-2) CONVG=.true.
        if (abs((wtaev - wtaev0)/wtaev) < 0.001_real64*accy) convg = .true.
        wtaev0 = wtaev
        if (.not.convg) cycle
        if (ldbpr(25) .and. .not.ldbpr(24) .and. mpi_rank==0) call print_radial_wavefunction (0)
        !IF (EOL) CALL matrix (dvdfirst)
        exit

    end do

    if (nit > nscf) then
        if (mpi_rank == 0) write (stderr_unit, *) " Maximum iterations in SCF Exceeded."
    end if

    do i = 31, 32 + max_mcp_file_multipole_rank
        close(i)                                ! The MCP coefficient files
    end do

    if (mpi_rank == 0) then
        !CLOSE (23)     ! The .rwf file
        close(mixing_output_unit)                               ! The .mix file
        !
        !   Complete the summary - moved from rscf92 for easier alloc/dalloc
        !
        call finish_summary
    end if
    !
    !   Deallocate storage
    !
    call dalloc (wt, "WT", "SCF")                       !Either getold or getald

    if (nec > 0) then
        call dalloc (iecc, "IECC", "SCF")
        call dalloc (ecv, "ECV", "SCF")
        !GG         CALL DALLOC (IQAR, 'IQAR', 'SCF')
        call dalloc (csf_subshell_occupations, "IQA", "SCF")
    end if

    if (nddim > 0) then
        call dalloc (da, "DA", "SCF")
        call dalloc (nda, "NDA", "SCF")
    end if

    if (nxdim > 0) then
        call dalloc (xa, "XA", "SCF")
        call dalloc (nxa, "NXA", "SCF")
    end if

    if (nydim > 0) then
        call dalloc (ya, "YA", "SCF")
        call dalloc (nya, "NYA", "SCF")
    end if

    if (extended_optimal_level) then
        call dalloc (shifted_eigenvalues, "EVAL", "SCF")
        call dalloc (csf_mixing_coefficients, "EvEC", "SCF")
        !GG         CALL DALLOC (IATJPO, 'IATJPO', 'SCF')
        !GG         CALL DALLOC (IASPAR, 'IASPAR', 'SCF')
        call dalloc (ncmaxblk, "NCMAXBLK", "SCF") ! getold.f
        call dalloc (eavblk, "EAVBLK", "SCF")     ! getold.f
        call dalloc (idxblk, "IDXBLK", "SCF")     ! Allocated in getold.f
        call dalloc (iccmin, "ICCMIN", "SCF")       ! Allocated in items.f<-getold.f
    end if
    !
    call dalloc (csf_offsets_by_block, "NCFPAST", "SCF")
    call dalloc (eigenvalue_offsets_by_block, "NCMINPAST", "SCF")
    call dalloc (mixing_offsets_by_block, "NEVECPAST", "SCF")


    return
end subroutine run_scf_iterations

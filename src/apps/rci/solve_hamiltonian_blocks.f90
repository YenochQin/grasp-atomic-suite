!***********************************************************************
!                                                                      *
subroutine solve_hamiltonian_blocks (unused_core_subshell_count, max_orbital_two_j)
    ! Module dependencies
    use calculation_file_units, only: mixing_output_unit, summary_output_unit
    !                                                                      *
    !   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
    !   and to diagonalise it.   The total angular momenta and parity of   *
    !   the ASF are found;  the eigenvectors are chosen so that the sign   *
    !   of the largest element is positive.                                *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
    !                        WGHTD5.                                       *
    !               [RCI92]: MANEIG, QED, SETHAM.                          *
    !                                                                      *
    !                                         Last revision: 28 Dec 1992   *
    !   Modified by Xinghong He               Last revision: 31 Dec 1997   *
    !   Block version Xinghong He             Last revision: 12 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use grasp_parameters,   only: max_orbitals
    use array_memory_management, only: alloc, dalloc
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use atomic_calculation_state, only: hartree_to_ev, hartree_to_wavenumber_cm, iccmin
    use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
        csf_mixing_coefficients, evecff, evecii
    !      USE foparm_C
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
        matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
        nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use angular_symmetry_labels, only: right_justified_j_labels => jlbr, parity_labels => jlbp
    use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
    use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
    use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
    use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use ci_restart_units, only: restart_matrix_unit
    use rci_prerun_settings, only: coeffcut1, coeffcut2, iprerun, ncsfpre
    use scf_previous_eigenvectors, only: evec1
    !USE blk_C
    use csf_block_cutoffs, only: iccut, iccutblk
    use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
    use coulomb_integral_table, only: frstco, indoei, ncoei, ndcoea, valoei
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
    use vacuum_polarization_integral_table, only: first_vacuum_polarization_lookup => frstvp, &
        & vacuum_polarization_integral_indices => indvpi, ndvpa, nvpi, vacuum_polarization_integral_values => valvpi
    use kinetic_energy_integral_table, only: frstki, indkei, ndkea, nkei, valkei
    use vinti_integral_table, only: first_vinti_lookup => frstvi, indtei, ndvin, nvinti, valtei
    use rci_matrix_file_position, only: nposition

    ! Procedure interfaces
    use prepare_correction_integral_tables_interface, only: prepare_correction_integral_tables

    ! Module dependencies
    use rci_csf_block_reader, only: load_csf_block

    ! Procedure interfaces
    use tabulate_self_energy_integrals_interface, only: tabulate_self_energy_integrals
    use assemble_ci_matrix_interface, only: assemble_ci_matrix
    use assemble_ci_matrix_from_restart_interface, only: assemble_ci_matrix_from_restart
    use write_hamiltonian_records_interface, only: write_hamiltonian_records

    ! Module dependencies
    use csf_quantum_queries, only: iq

    ! Procedure interfaces
    use diagonalize_hamiltonian_block_interface, only: diagonalize_hamiltonian_block
    use write_level_energies_interface, only: write_level_energies
    use read_ci_level_selection_interface, only: read_ci_level_selection
    use build_qed_self_energy_interface, only: build_qed_self_energy
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: unused_core_subshell_count
    integer,  intent(in):: max_orbital_two_j

    ! Local declarations
    real(real64) :: slfint(max_orbitals)
    character(len=8) :: number_text
    real(real64) :: atwinv, elsto, eau, ecm, eev, elemnt
    real(real64), pointer :: slf_en(:), ucf(:), etot(:)
    integer(int64) :: nelmnt_a
    integer :: mpi_rank, mpi_process_count,  iiatjpo, iiaspar
    integer :: i, j, irestart, ncminpas, symmetry_block_index, csf_index, ip, mode
    !-----------------------------------------------
    !
    mpi_rank = 0
    mpi_process_count = 1
    !***************************************************************
    !     ...Common to all blocks - place here to save CPU time
    call prepare_correction_integral_tables (max_orbital_two_j, atwinv)

    !***************************************************************
    !      Loop over blocks
    !***************************************************************
    ncminpas = 0

    do symmetry_block_index = 1, symmetry_block_count
        csf_count    =   csf_counts_by_block(symmetry_block_index)
        nvec   =   eigenvector_counts_by_block(symmetry_block_index)
        nvecmx = ncmaxblk(symmetry_block_index)
        iccut(1)  = iccutblk(symmetry_block_index)
        call alloc (slf_en,csf_count,"SLF_EN", "MATRIX")
        call alloc (ucf, orbital_count,"UCF", "MATRIX")
        do csf_index=1,csf_count
            slf_en(csf_index) = 0.0
        end do

        nposition = 7 + orbital_count + orbital_count  ! File position of the previous block
        ! in the .res file
        do i = 1, symmetry_block_index - 1
            j = (csf_counts_by_block(i) - mpi_rank - 1 + mpi_process_count) / mpi_process_count
            if (csf_counts_by_block(i) < mpi_process_count) j = csf_counts_by_block(i) / (mpi_rank+1)
            nposition = nposition + j + 1
        end do

        !.. SETHAM does not handle this extrem case
        if (mpi_process_count > csf_count)                                          &
            stop "matrix: too many nodes"

        !        ...Obtain ivec() from iccmin()
        if (nvec > 0) then
            call alloc (ivec, nvec, "IVEC", "MATRIX")
            do i = 1, nvec
                ivec(i) = iccmin(i+ncminpas)
            end do
            ncminpas = ncminpas + nvec
        end if

        !        ...These 3 were allocated in lodcsh2 and deallocated at the end
        !        ... of this routine and in the setham. In this block version,
        !        ... both allocation and deallocation are placed here.
        call alloc (csf_subshell_occupations, max_orbitals, csf_count, "IQA", "MATRIX")
        call alloc (csf_subshell_quantum_codes, max_orbitals, 3, csf_count, "JQSA", "MATRIX")
        call alloc (csf_coupling_quantum_codes, max_orbitals, csf_count, "JCUPA", "MATRIX")

        !      ...Load CSF list of the current block
        block
            ! Local declarations
            integer :: input_status
            call load_csf_block(symmetry_block_index,input_status)
            if (input_status /= 0) error stop "RCI descriptor block failed"
        end block
        ! zou

        if (lse) then
            print *, "Entering QED ..."

            write (summary_output_unit,*)
            write (summary_output_unit,*) " Self Energy Corrections: "
            write (summary_output_unit,*)
            call tabulate_self_energy_integrals (slfint)
            do csf_index = 1, csf_count
                elemnt = 0.0_real64
                do i = 1,orbital_count
                    !              print *,ic, i,iq(i,ic),slfint(i)
                    elemnt = elemnt+iq (i,csf_index)*slfint(i)
                end do
                !              print *,ic, elemnt
                slf_en(csf_index) = elemnt
            end do

            write (summary_output_unit,*)
            write (summary_output_unit,*) "Self-energy corrections estimated"           &
                //" --- these will influence the data"
            write (summary_output_unit,*) " in the RCI92 MIXing coefficients File."
        end if

        ! zou
        if (nvec <= 0) then
            write (mixing_output_unit) symmetry_block_index, csf_count, nvec, 999, 999
            !           ...Generate H matrix of the current block and write to disk
            !           ...eav, nelmnt are also obtained from genmat
            call assemble_ci_matrix (atwinv, symmetry_block_index, mpi_rank, mpi_process_count, elsto, irestart,&
                slf_en)
            ! This call is optional, meaning the matrix of this block
            ! does not need to be computed. But don't comment it out
            ! since other programs may assume thet existence of them.
            call assemble_ci_matrix_from_restart (irestart, nelmnt_a, elsto)
        else
            call assemble_ci_matrix (atwinv, symmetry_block_index, mpi_rank, mpi_process_count, elsto, irestart,   &
                slf_en)
            call assemble_ci_matrix_from_restart (irestart, nelmnt_a, elsto)
            !
            !   Allocate and deallocate memory for the mixing coefficients
            !   from the prerun
            !
            if (iprerun==1) call alloc (evec1,csf_count*nvec, "EVEC1", "MATRIX")
            if (iprerun==2) call dalloc (evec1, "EVEC1", "MATRIX")

            call diagonalize_hamiltonian_block (iiatjpo, iiaspar)
            !
            !   Write out eigenvalues (ENGOUT), dominant components of the
            !   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
            !   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
            !   EAV and ELSTO are added back to energy here
            !
            ! ELSTO has never been in Hamiltonian matrix, yet it was
            ! added to EAV which was later substracted from H. Thus at
            ! this point, EAV is correct (it has ELSTO added), EVAL()
            ! need ELSTO and the correct EAV.
            if (csf_count > 1) then
                do i = 1, nvec
                    shifted_eigenvalues(i) = shifted_eigenvalues(i) + elsto
                end do
            end if

            call write_level_energies (hamiltonian_energy_shift,shifted_eigenvalues,iiatjpo,iiaspar,ivec,nvec,3)
            call read_ci_level_selection (iiatjpo, iiaspar)

            !     ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
            !     ...MIXing coefficients File; close the file; print a report
            write (mixing_output_unit) symmetry_block_index, csf_count, nvec, iiatjpo, iiaspar
            write (mixing_output_unit) (ivec(i), i = 1, nvec)
            write (mixing_output_unit) hamiltonian_energy_shift,(shifted_eigenvalues(i),i = 1,nvec)
            write (mixing_output_unit) ((csf_mixing_coefficients(i+(j-1)*csf_count),i = 1,csf_count),j = 1,nvec)
            !      print '(5i5)',jblock, ncf, nvec, iiatjpo, iiaspar
            !      print '(3I20)',(ivec(i), i = 1, nvec)
            !      print '(3f20.15)',EAV,(EVAL(I),I = 1,NVEC)
            !      print '(3f20.15)',((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
            print *, "RCI92 MIXing coefficients File generated."
            !pdbg
            !   Save the mixing coefficients from the prerun
            !
            if (iprerun == 1) then
                do j = 1, nvec
                    do i = 1, csf_count
                        evec1(i+(j-1)*csf_count) = csf_mixing_coefficients(i+(j-1)*csf_count)
                    end do
                end do
            end if
            !
            !   Estimate diagonal self-energy corrections; the stored
            !   eigenvalues and eigenvectors are not modified by these
            !   estimates
            !
            if (.not.lse) then
                print *, "Entering QED ..."
                call alloc (etot,nvec,"ETOT", "MATRIX")

                write (summary_output_unit,*)
                write (summary_output_unit,*) " Self Energy Corrections: "
                write (summary_output_unit,*)
                write (summary_output_unit,'(" Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV")')
                write (summary_output_unit,*)
                do j = 1, nvec
                    call build_qed_self_energy (j,slfint,ucf)
                    elemnt = 0.0_real64
                    csf_index = ivec(j)
                    do i = 1,orbital_count
                        elemnt = elemnt+ucf(i)*slfint(i)
                        !                 ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
                    end do
                    etot(j) = shifted_eigenvalues(j)+elemnt
                    !
                    eau = elemnt
                    ecm = eau*hartree_to_wavenumber_cm
                    eev = eau*hartree_to_ev
                    ip = (iiaspar+3)/2
                    write (summary_output_unit,'(1i3,2x,2a4,1p,3d22.14)') j,right_justified_j_labels(iiatjpo),parity_labels( &
                        & ip),eau,ecm,eev
                    !
                end do

                write (summary_output_unit,*)
                write (summary_output_unit,*) "Self-energy corrections estimated"           &
                    //" --- these do not influence the data"
                write (summary_output_unit,*) " in the RCI92 MIXing coefficients File."
                call write_level_energies (hamiltonian_energy_shift,etot,iiatjpo,iiaspar,ivec,nvec,mode)
                ! zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
                call dalloc (etot, "ETOT", "MATRIX")
            end if

            !        ...Locals
            call dalloc (ivec, "IVEC", "MATRIX")
            !        ...Allocated in maneig
            call dalloc (shifted_eigenvalues, "EVAL", "MATRIX")
            call dalloc (csf_mixing_coefficients, "EVEC", "MATRIX")

        end if

        !        ...Locals
        call dalloc (csf_subshell_occupations, "IQA", "MATRIX")
        call dalloc (csf_subshell_quantum_codes, "JQSA", "MATRIX")
        call dalloc (csf_coupling_quantum_codes, "JCUPA", "MATRIX")
        call dalloc (slf_en, "SLF_EN", "MATRIX")
        call dalloc (ucf, "UCF", "MATRIX")

    end do
    !
    !   Close the restart files; nothing will be added to them now
    !
    close (restart_matrix_unit)
    !     ...Clean up
    call dalloc (csf_counts_by_block, "NCFBLK", "MATRIX")
    call dalloc (eigenvector_counts_by_block, "NEVBLK", "MATRIX")
    call dalloc (ncmaxblk, "NCMAXBLK", "MATRIX")
    call dalloc (iccutblk, "ICUTTBLK", "MATRIX")
    call dalloc (iccmin, "ICCMIN", "MATRIX") ! allocated in items as pnccmin

    call dalloc (valteirk, "VALTEIRK", "MATRIX") ! allocated in genintrk
    call dalloc (indteirk, "INDTEIRK", "MATRIX") ! allocated in genintrk
    !
    !   Deallocate storage for the integral lists from the
    !   Dirac-Coulomb operator; the storage was allocated
    !   in IABINT and RKINTC
    !
    if (ncoei > 0) then
        call dalloc (indoei, "INDOEI", "MATRIX")
        call dalloc (valoei, "VALOEI", "MATRIX")
    end if
    !
    !   Deallocate storage for the integral lists from the
    !   transverse photon interaction operator; this storage
    !   was allocated in BRINT1, brint2,...brint6
    !
    if (ltrans) then
        if (ntpi(1) > 0) then
            call dalloc (indtp1, "INDTP1", "MATRIX")
            call dalloc (valtp1, "VALTP1", "MATRIX")
        end if
        if (ntpi(2) > 0) then
            call dalloc (indtp2, "INDTP2", "MATRIX")
            call dalloc (valtp2, "VALTP2", "MATRIX")
        end if
        if (ntpi(3) > 0) then
            call dalloc (indtp3, "INDTP3", "MATRIX")
            call dalloc (valtp3, "VALTP3", "MATRIX")
        end if
        if (ntpi(4) > 0) then
            call dalloc (indtp4, "INDTP4", "MATRIX")
            call dalloc (valtp4, "VALTP4", "MATRIX")
        end if
        if (ntpi(5) > 0) then
            call dalloc (indtp5, "INDTP5", "MATRIX")
            call dalloc (valtp5, "VALTP5", "MATRIX")
        end if
        if (ntpi(6) > 0) then
            call dalloc (indtp6, "INDTP6", "MATRIX")
            call dalloc (valtp6, "VALTP6", "MATRIX")
        end if
    end if
    !
    !   Deallocate storage for the nuclear motional energy integral
    !   lists; this was allocated in KEINT and VINT
    !
    if (lnms) then
        if (nkei > 0) then
            call dalloc (indkei, "INDKEI", "MATRIX")
            call dalloc (valkei, "VALKEI", "MATRIX")
        end if
    end if
    if (lsms) then
        if (nvinti > 0) then
            call dalloc (indtei, "INDTEI", "MATRIX")
            call dalloc (valtei, "VALTIE", "MATRIX")
        end if
    end if
    !
    !   Deallocate storage for the vacuum-polarisation integral list;
    !   this was allocated in VPINT
    !
    if (lvp) then
        if (nvpi > 0) then
            call dalloc (vacuum_polarization_integral_indices, "INDVPI", "MATRIX")
            call dalloc (vacuum_polarization_integral_values, "VALVPI", "MATRIX")
        end if
    end if

    call dalloc (radial_large_components, "PF", "MATRIX") ! lodrwf or lodres
    call dalloc (radial_small_components, "QF", "MATRIX") ! lodrwf or lodres

    return

end subroutine solve_hamiltonian_blocks

module solve_hamiltonian_blocks_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: solve_hamiltonian_blocks
contains

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
        use hamiltonian_correction_flags, only: include_normal_mass_shift, estimate_self_energy, include_specific_mass_shift, &
            include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state, only: hartree_to_ev, hartree_to_wavenumber_cm, requested_level_indices
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        !      USE foparm_C
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, &
            highest_selected_level_by_block, symmetry_block_count
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
        use selected_level_indices, only: active_level_indices, active_level_count, highest_active_level_index
        use csf_coupling_data, only: csf_angular_codes, csf_subshell_quantum_codes, csf_angular_code_count
        use radial_wavefunctions, only: max_orbitals, radial_large_components, radial_small_components
        use ci_restart_units, only: restart_matrix_unit
        use rci_prerun_settings, only: prerun_stage
        use scf_previous_eigenvectors, only: previous_csf_mixing_coefficients
        !USE blk_C
        use csf_block_cutoffs, only: active_csf_cutoffs, csf_cutoff_by_block
        use distributed_slater_integrals, only: local_slater_integral_indices, local_slater_integral_values
        use coulomb_integral_table, only: coulomb_integral_indices, coulomb_integral_count, coulomb_integral_values
        use breit_integral_tables, only: breit_integral_counts, breit_integral_indices_1, breit_integral_indices_2, &
            breit_integral_indices_3, breit_integral_indices_4, breit_integral_indices_5, breit_integral_indices_6, &
            breit_integral_values_1, breit_integral_values_2, breit_integral_values_3, breit_integral_values_4, &
            breit_integral_values_5, breit_integral_values_6
        use vacuum_polarization_integral_table, only: vacuum_polarization_integral_indices, &
            vacuum_polarization_integral_count, vacuum_polarization_integral_values
        use kinetic_energy_integral_table, only: kinetic_integral_indices, kinetic_integral_count, kinetic_integral_values
        use vinti_integral_table, only: vinti_integral_indices, vinti_integral_count, vinti_integral_values
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
        use csf_quantum_queries, only: csf_occupation

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
        real(real64) :: orbital_self_energy_integrals(max_orbitals)
        real(real64) :: inverse_nuclear_mass, closed_shell_breit_shift, energy_hartree, energy_wavenumber, energy_ev, &
            & self_energy_shift
        real(real64), pointer :: csf_self_energies(:), orbital_occupation_weights(:), corrected_level_energies(:)
        integer(int64) :: nelmnt_a
        integer :: mpi_rank, mpi_process_count,  block_two_j_plus_one, block_parity_sign
        integer :: previous_block_index, local_csf_record_count, level_index, coefficient_csf_index, orbital_index
        integer :: restart_header_status, selected_level_offset, symmetry_block_index, csf_index, parity_label_index, &
            & energy_difference_mode
        !-----------------------------------------------
        !
        mpi_rank = 0
        mpi_process_count = 1
        !***************************************************************
        !     ...Common to all blocks - place here to save CPU time
        call prepare_correction_integral_tables (max_orbital_two_j, inverse_nuclear_mass)

        !***************************************************************
        !      Loop over blocks
        !***************************************************************
        selected_level_offset = 0

        do symmetry_block_index = 1, symmetry_block_count
            csf_count    =   csf_counts_by_block(symmetry_block_index)
            active_level_count   =   eigenvector_counts_by_block(symmetry_block_index)
            highest_active_level_index = highest_selected_level_by_block(symmetry_block_index)
            active_csf_cutoffs(1)  = csf_cutoff_by_block(symmetry_block_index)
            call alloc (csf_self_energies,csf_count,"csf_self_energies", "solve_hamiltonian_blocks")
            call alloc (orbital_occupation_weights, orbital_count,"orbital_occupation_weights", "solve_hamiltonian_blocks")
            do csf_index=1,csf_count
                csf_self_energies(csf_index) = 0.0
            end do

            nposition = 7 + orbital_count + orbital_count  ! File position of the previous block
            ! in the .res file
            do previous_block_index = 1, symmetry_block_index - 1
                local_csf_record_count = (csf_counts_by_block(previous_block_index) - mpi_rank - 1 + mpi_process_count) / &
                    & mpi_process_count
                if (csf_counts_by_block(previous_block_index) < mpi_process_count) local_csf_record_count = &
                    & csf_counts_by_block(previous_block_index) / (mpi_rank+1)
                nposition = nposition + local_csf_record_count + 1
            end do

            !.. SETHAM does not handle this extrem case
            if (mpi_process_count > csf_count)                                          &
                stop "solve_hamiltonian_blocks: too many nodes"

            !        ...Obtain ivec() from iccmin()
            if (active_level_count > 0) then
                call alloc (active_level_indices, active_level_count, "active_level_indices", "solve_hamiltonian_blocks")
                do level_index = 1, active_level_count
                    active_level_indices(level_index) = requested_level_indices(level_index+selected_level_offset)
                end do
                selected_level_offset = selected_level_offset + active_level_count
            end if

            !        ...These 3 were allocated in lodcsh2 and deallocated at the end
            !        ... of this routine and in the setham. In this block version,
            !        ... both allocation and deallocation are placed here.
            call alloc (csf_subshell_occupations, max_orbitals, csf_count, "csf_subshell_occupations", "solve_hamiltonian_blocks")
            call alloc (csf_subshell_quantum_codes, max_orbitals, 3, csf_count, "csf_subshell_quantum_codes", &
                "solve_hamiltonian_blocks")
            call alloc (csf_angular_codes, csf_angular_code_count, csf_count, "csf_angular_codes", &
                "solve_hamiltonian_blocks")

            !      ...Load CSF list of the current block
            block
                ! Local declarations
                integer :: input_status
                call load_csf_block(symmetry_block_index,input_status)
                if (input_status /= 0) error stop "RCI descriptor block failed"
            end block
            ! zou

            if (estimate_self_energy) then
                print *, "Entering QED ..."

                write (summary_output_unit,*)
                write (summary_output_unit,*) " Self Energy Corrections: "
                write (summary_output_unit,*)
                call tabulate_self_energy_integrals (orbital_self_energy_integrals)
                do csf_index = 1, csf_count
                    self_energy_shift = 0.0_real64
                    do orbital_index = 1,orbital_count
                        !              print *,ic, orbital_index,iq(orbital_index,ic),slfint(orbital_index)
                        self_energy_shift = self_energy_shift+csf_occupation &
                            & (orbital_index,csf_index)*orbital_self_energy_integrals(orbital_index)
                    end do
                    !              print *,ic, elemnt
                    csf_self_energies(csf_index) = self_energy_shift
                end do

                write (summary_output_unit,*)
                write (summary_output_unit,*) "Self-energy corrections estimated"           &
                    //" --- these will influence the data"
                write (summary_output_unit,*) " in the RCI92 MIXing coefficients File."
            end if

            ! zou
            if (active_level_count <= 0) then
                write (mixing_output_unit) symmetry_block_index, csf_count, active_level_count, 999, 999
                !           ...Generate H matrix of the current block and write to disk
                !           ...eav, nelmnt are also obtained from genmat
                call assemble_ci_matrix (inverse_nuclear_mass, symmetry_block_index, mpi_rank, mpi_process_count, &
                    & closed_shell_breit_shift, restart_header_status,&
                    csf_self_energies)
                ! This call is optional, meaning the matrix of this block
                ! does not need to be computed. But don't comment it out
                ! since other programs may assume thet existence of them.
                call assemble_ci_matrix_from_restart (restart_header_status, nelmnt_a, closed_shell_breit_shift)
            else
                call assemble_ci_matrix (inverse_nuclear_mass, symmetry_block_index, mpi_rank, mpi_process_count, &
                    & closed_shell_breit_shift, restart_header_status,   &
                    csf_self_energies)
                call assemble_ci_matrix_from_restart (restart_header_status, nelmnt_a, closed_shell_breit_shift)
                !
                !   Allocate and deallocate memory for the mixing coefficients
                !   from the prerun
                !
                if (prerun_stage==1) call alloc (previous_csf_mixing_coefficients,csf_count*active_level_count, &
                    "previous_csf_mixing_coefficients", "solve_hamiltonian_blocks")
                if (prerun_stage==2) call dalloc (previous_csf_mixing_coefficients, "previous_csf_mixing_coefficients", &
                    "solve_hamiltonian_blocks")

                call diagonalize_hamiltonian_block (block_two_j_plus_one, block_parity_sign)
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
                    do level_index = 1, active_level_count
                        shifted_eigenvalues(level_index) = shifted_eigenvalues(level_index) + closed_shell_breit_shift
                    end do
                end if

                call write_level_energies &
                    & (hamiltonian_energy_shift,shifted_eigenvalues,block_two_j_plus_one,block_parity_sign,active_level_indices, &
                    active_level_count,3)
                call read_ci_level_selection (block_two_j_plus_one, block_parity_sign)

                !     ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
                !     ...MIXing coefficients File; close the file; print a report
                write (mixing_output_unit) symmetry_block_index, csf_count, active_level_count, block_two_j_plus_one, &
                    & block_parity_sign
                write (mixing_output_unit) (active_level_indices(level_index), level_index = 1, active_level_count)
                write (mixing_output_unit) hamiltonian_energy_shift,(shifted_eigenvalues(level_index),level_index = &
                    & 1,active_level_count)
                write (mixing_output_unit) &
                    & ((csf_mixing_coefficients(coefficient_csf_index+(level_index-1)*csf_count),coefficient_csf_index = &
                    & 1,csf_count),level_index = 1,active_level_count)
                !      print '(5i5)',jblock, ncf, nvec, iiatjpo, iiaspar
                !      print '(3I20)',(ivec(i), i = 1, nvec)
                !      print '(3f20.15)',EAV,(EVAL(I),I = 1,NVEC)
                !      print '(3f20.15)',((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)
                print *, "RCI92 MIXing coefficients File generated."
                !pdbg
                !   Save the mixing coefficients from the prerun
                !
                if (prerun_stage == 1) then
                    do level_index = 1, active_level_count
                        do coefficient_csf_index = 1, csf_count
                            previous_csf_mixing_coefficients(coefficient_csf_index+(level_index-1)*csf_count) = &
                                & csf_mixing_coefficients(coefficient_csf_index+(level_index-1)*csf_count)
                        end do
                    end do
                end if
                !
                !   Estimate diagonal self-energy corrections; the stored
                !   eigenvalues and eigenvectors are not modified by these
                !   estimates
                !
                if (.not.estimate_self_energy) then
                    print *, "Entering QED ..."
                    call alloc (corrected_level_energies,active_level_count,"corrected_level_energies", "solve_hamiltonian_blocks")

                    write (summary_output_unit,*)
                    write (summary_output_unit,*) " Self Energy Corrections: "
                    write (summary_output_unit,*)
                    write (summary_output_unit,'(" Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV")')
                    write (summary_output_unit,*)
                    do level_index = 1, active_level_count
                        call build_qed_self_energy (level_index,orbital_self_energy_integrals,orbital_occupation_weights)
                        self_energy_shift = 0.0_real64
                        csf_index = active_level_indices(level_index)
                        do orbital_index = 1,orbital_count
                            self_energy_shift = &
                                & self_energy_shift + orbital_occupation_weights(orbital_index) &
                                & *orbital_self_energy_integrals(orbital_index)
                            !                 ELEMNT = ELEMNT+IQ (I,IC)*SLFINT(I)
                        end do
                        corrected_level_energies(level_index) = shifted_eigenvalues(level_index)+self_energy_shift
                        !
                        energy_hartree = self_energy_shift
                        energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
                        energy_ev = energy_hartree*hartree_to_ev
                        parity_label_index = (block_parity_sign+3)/2
                        write (summary_output_unit,'(1i3,2x,2a,1p,3d22.14)') &
                            & level_index,format_j_label(block_two_j_plus_one),parity_labels( &
                            & parity_label_index),energy_hartree,energy_wavenumber,energy_ev
                        !
                    end do

                    write (summary_output_unit,*)
                    write (summary_output_unit,*) "Self-energy corrections estimated"           &
                        //" --- these do not influence the data"
                    write (summary_output_unit,*) " in the RCI92 MIXing coefficients File."
                    call write_level_energies &
                        & (hamiltonian_energy_shift, corrected_level_energies, block_two_j_plus_one, &
                            block_parity_sign, active_level_indices, &
                        active_level_count,energy_difference_mode)
                    ! zou       CALL ENGOUT (EAV+elsto,ETOT,IiATJPO,iIASPAR,IVEC,NVEC,MODE)
                    call dalloc (corrected_level_energies, "corrected_level_energies", "solve_hamiltonian_blocks")
                end if

                !        ...Locals
                call dalloc (active_level_indices, "active_level_indices", "solve_hamiltonian_blocks")
                !        ...Allocated in maneig
                call dalloc (shifted_eigenvalues, "shifted_eigenvalues", "solve_hamiltonian_blocks")
                call dalloc (csf_mixing_coefficients, "csf_mixing_coefficients", "solve_hamiltonian_blocks")

            end if

            !        ...Locals
            call dalloc (csf_subshell_occupations, "csf_subshell_occupations", "solve_hamiltonian_blocks")
            call dalloc (csf_subshell_quantum_codes, "csf_subshell_quantum_codes", "solve_hamiltonian_blocks")
            call dalloc (csf_angular_codes, "csf_angular_codes", "solve_hamiltonian_blocks")
            call dalloc (csf_self_energies, "csf_self_energies", "solve_hamiltonian_blocks")
            call dalloc (orbital_occupation_weights, "orbital_occupation_weights", "solve_hamiltonian_blocks")

        end do
        !
        !   Close the restart files; nothing will be added to them now
        !
        close (restart_matrix_unit)
        !     ...Clean up
        call dalloc (csf_counts_by_block, "csf_counts_by_block", "solve_hamiltonian_blocks")
        call dalloc (eigenvector_counts_by_block, "eigenvector_counts_by_block", "solve_hamiltonian_blocks")
        call dalloc (highest_selected_level_by_block, "highest_selected_level_by_block", "solve_hamiltonian_blocks")
        call dalloc (csf_cutoff_by_block, "csf_cutoff_by_block", "solve_hamiltonian_blocks")
        call dalloc (requested_level_indices, "requested_level_indices", "solve_hamiltonian_blocks") ! allocated in items as pnccmin

        call dalloc (local_slater_integral_values, "local_slater_integral_values", "solve_hamiltonian_blocks") ! allocated in genintrk
        call dalloc (local_slater_integral_indices, "local_slater_integral_indices", "solve_hamiltonian_blocks") ! allocated in genintrk
        !
        !   Deallocate storage for the integral lists from the
        !   Dirac-Coulomb operator; the storage was allocated
        !   in IABINT and RKINTC
        !
        if (coulomb_integral_count > 0) then
            call dalloc (coulomb_integral_indices, "coulomb_integral_indices", "solve_hamiltonian_blocks")
            call dalloc (coulomb_integral_values, "coulomb_integral_values", "solve_hamiltonian_blocks")
        end if
        !
        !   Deallocate storage for the integral lists from the
        !   transverse photon interaction operator; this storage
        !   was allocated in BRINT1, brint2,...brint6
        !
        if (include_transverse_interaction) then
            if (breit_integral_counts(1) > 0) then
                call dalloc (breit_integral_indices_1, "breit_integral_indices_1", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_1, "breit_integral_values_1", "solve_hamiltonian_blocks")
            end if
            if (breit_integral_counts(2) > 0) then
                call dalloc (breit_integral_indices_2, "breit_integral_indices_2", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_2, "breit_integral_values_2", "solve_hamiltonian_blocks")
            end if
            if (breit_integral_counts(3) > 0) then
                call dalloc (breit_integral_indices_3, "breit_integral_indices_3", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_3, "breit_integral_values_3", "solve_hamiltonian_blocks")
            end if
            if (breit_integral_counts(4) > 0) then
                call dalloc (breit_integral_indices_4, "breit_integral_indices_4", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_4, "breit_integral_values_4", "solve_hamiltonian_blocks")
            end if
            if (breit_integral_counts(5) > 0) then
                call dalloc (breit_integral_indices_5, "breit_integral_indices_5", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_5, "breit_integral_values_5", "solve_hamiltonian_blocks")
            end if
            if (breit_integral_counts(6) > 0) then
                call dalloc (breit_integral_indices_6, "breit_integral_indices_6", "solve_hamiltonian_blocks")
                call dalloc (breit_integral_values_6, "breit_integral_values_6", "solve_hamiltonian_blocks")
            end if
        end if
        !
        !   Deallocate storage for the nuclear motional energy integral
        !   lists; this was allocated in KEINT and VINT
        !
        if (include_normal_mass_shift) then
            if (kinetic_integral_count > 0) then
                call dalloc (kinetic_integral_indices, "kinetic_integral_indices", "solve_hamiltonian_blocks")
                call dalloc (kinetic_integral_values, "kinetic_integral_values", "solve_hamiltonian_blocks")
            end if
        end if
        if (include_specific_mass_shift) then
            if (vinti_integral_count > 0) then
                call dalloc (vinti_integral_indices, "vinti_integral_indices", "solve_hamiltonian_blocks")
                call dalloc (vinti_integral_values, "vinti_integral_values", "solve_hamiltonian_blocks")
            end if
        end if
        !
        !   Deallocate storage for the vacuum-polarisation integral list;
        !   this was allocated in VPINT
        !
        if (include_vacuum_polarization) then
            if (vacuum_polarization_integral_count > 0) then
                call dalloc (vacuum_polarization_integral_indices, "vacuum_polarization_integral_indices", &
                    "solve_hamiltonian_blocks")
                call dalloc (vacuum_polarization_integral_values, "vacuum_polarization_integral_values", "solve_hamiltonian_blocks")
            end if
        end if

        call dalloc (radial_large_components, "radial_large_components", "solve_hamiltonian_blocks") ! lodrwf or lodres
        call dalloc (radial_small_components, "radial_small_components", "solve_hamiltonian_blocks") ! lodrwf or lodres

        return

    end subroutine solve_hamiltonian_blocks
end module solve_hamiltonian_blocks_interface

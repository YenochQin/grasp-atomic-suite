! Restore the saved CI matrix blocks for diagonalization.
module assemble_ci_matrix_from_restart_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: assemble_ci_matrix_from_restart
contains

    !   Get eav and do writings to the summary file .csum
    !   The mpi version (genmat2mpi) also gets nelmnt_a and elsto
    !
    ! Xinghong He 98-06-15
    subroutine assemble_ci_matrix_from_restart(restart_header_status, nelmnt_a, closed_shell_breit_shift)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        use, intrinsic :: iso_fortran_env, only: real64, int64
        use sparse_hamiltonian_storage, only: block_matrix_element_cutoff, block_assembly_csf_count, &
            block_coulomb_one_body_coefficient_count, block_coulomb_one_body_integral_count, &
            block_breit_core_coefficient_count, block_coulomb_two_body_coefficient_count, block_slater_integral_count, &
            block_stored_element_count, block_kinetic_integral_count, block_breit_coefficient_count, &
            block_breit_integral_counts, block_vinti_integral_count, block_vacuum_polarization_integral_count
        use hamiltonian_correction_flags, only: include_normal_mass_shift, include_specific_mass_shift, &
            include_transverse_interaction, include_vacuum_polarization
        use eigenstate_expansions, only: hamiltonian_energy_shift
        use mpi_state, only: mpi_allreduce, mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, &
            mpi_integer8, mpi_max, mpi_rank, mpi_reduce, mpi_sum
        implicit none (type, external)

        ! Arguments
        integer, intent(in)        :: restart_header_status
        integer(int64), intent(out) :: nelmnt_a
        real(real64), intent(in)   :: closed_shell_breit_shift

        ! Local declarations
        integer :: ntpi_a(6)
        integer               :: ncoei_a, nctei_a, nctec_a, nmcbp_a,    &
            ncore_a, nvpi_a, nkei_a, nvinti_a,     &
            ncoec_a
        real(real64)          :: eav_a
        ! Transfer parameters around and print out on node-0.
        !  .Temporary variables (ending with "_a") are used to store the
        !   sums over nodes;
        !  .Outputs to stream 24 are done on node-0;
        ! The following parameters are accumulated in setham which will not
        ! contain the correct values in restart mode. And thus is skipped.
        if (restart_header_status /= 0) then                    ! non-restart mode

            call mpi_reduce (block_coulomb_one_body_integral_count, ncoei_a, 1, mpi_integer, mpi_sum, 0,  &
                mpi_comm_world, mpi_error_code)
            call mpi_reduce (block_coulomb_one_body_coefficient_count, ncoec_a, 1, mpi_integer, mpi_sum, 0,  &
                mpi_comm_world, mpi_error_code)
            call mpi_reduce (block_slater_integral_count, nctei_a, 1, mpi_integer, mpi_sum, 0,  &
                mpi_comm_world, mpi_error_code)
            call mpi_reduce (block_coulomb_two_body_coefficient_count, nctec_a, 1, mpi_integer, mpi_sum, 0,  &
                mpi_comm_world, mpi_error_code)

            if (include_transverse_interaction) then
                call mpi_reduce (block_breit_integral_counts, ntpi_a, 6, mpi_integer, mpi_sum, 0, &
                    mpi_comm_world, mpi_error_code)
                call mpi_reduce (block_breit_coefficient_count, nmcbp_a, 1, mpi_integer, mpi_sum,0,&
                    mpi_comm_world, mpi_error_code)
                call mpi_reduce (block_breit_core_coefficient_count, ncore_a, 1, mpi_integer,mpi_sum,0, &
                    mpi_comm_world, mpi_error_code)
            end if

            if (include_vacuum_polarization)                                                         &
                call mpi_reduce (block_vacuum_polarization_integral_count, nvpi_a, 1, mpi_integer, mpi_sum, 0, &
                mpi_comm_world, mpi_error_code)

            if (include_normal_mass_shift)                                                        &
                call mpi_reduce (block_kinetic_integral_count, nkei_a, 1, mpi_integer, mpi_sum, 0, &
                mpi_comm_world, mpi_error_code)

            if (include_specific_mass_shift)                                                        &
                call mpi_reduce (block_vinti_integral_count, nvinti_a, 1, mpi_integer, mpi_sum,&
                0, mpi_comm_world, mpi_error_code)

            if (mpi_rank == 0) then
                write (summary_output_unit,'("CUTOFF set to ",1p,d17.10)') block_matrix_element_cutoff
                write (summary_output_unit,'("Dirac-Coulomb one-e radial integrals:",1i8)') ncoei_a
                write (summary_output_unit,'("One-e angular integrals that exceed CUTOFF: ",1i8)') ncoec_a
                write (summary_output_unit,'("Coulomb two-e radial integrals: ",1i8)') nctei_a
                write (summary_output_unit,'("Two-e angular integrals that exceed CUTOFF: ",1i11)') nctec_a
                if (include_transverse_interaction) then
                    write (summary_output_unit,'("Transverse two-e radial integrals: "/,6i8)') ntpi_a
                    write (summary_output_unit,'("MCBP coefficients that exceed CUTOFF: ",1i16)') nmcbp_a
                    write (summary_output_unit,'("Core coefficients that exceed CUTOFF: ",1i8)') ncore_a
                end if
                if (include_vacuum_polarization) write (summary_output_unit,'("Vacuum polarisation integrals: ",1i8)') nvpi_a
                if (include_normal_mass_shift) write (summary_output_unit,'("Kinetic energy integrals: ",1i8)') nkei_a
                if (include_specific_mass_shift) write (summary_output_unit,'("Vinti integrals: ",1i8)') nvinti_a
            end if
        else
            if (mpi_rank == 0) then
                write (summary_output_unit,*) "Restart mode --- no report on radial integrals"
            end if
        end if                                      !(irestart .NE. 0) ! non-restart mode

        ! ELSTO, EAV are not only for print-out, but also used later.
        ! density of the Hamiltonian matrix is only for print-out.
        ! At this place Hamiltonian matrix elements (EMT on each node) do
        !   _not_ contain ELSTO. ELSTO will be added to the total energy
        !  later with EAV.

        call mpi_allreduce (block_stored_element_count, nelmnt_a, 1, mpi_integer8,  &
            ! AK:      MPI_INTEGER8 required to match the sizes of NELMNT_a and NELMNTtmp
            ! CFF &    MPI_SUM, MPI_COMM_WORLD, ierr) ! want MAX not SUM
            &    mpi_max, mpi_comm_world, mpi_error_code)
        ! NELMNT_a decides between sparse and dense
        ! see comments in maneig.f90

        call mpi_bcast (closed_shell_breit_shift, 1, mpi_double_precision, 0,         &
            mpi_comm_world, mpi_error_code)
        call mpi_allreduce (hamiltonian_energy_shift, eav_a, 1, mpi_double_precision,   &
            mpi_sum, mpi_comm_world, mpi_error_code)
        hamiltonian_energy_shift = eav_a / dble (block_assembly_csf_count) + closed_shell_breit_shift

        if (mpi_rank == 0) then
            write (summary_output_unit,'("Elements that exceed CUTOFF in the lower"," triangle of the H matrix: ",1i11)') nelmnt_a
            write (summary_output_unit, *)
            write (summary_output_unit, '("Average energy = ",1p,d19.12," Hartrees.")') hamiltonian_energy_shift
        end if

        ! cjb 1I8 -> 1I16

        return
    end subroutine assemble_ci_matrix_from_restart
end module assemble_ci_matrix_from_restart_interface

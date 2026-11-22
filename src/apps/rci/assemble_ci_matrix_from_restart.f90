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
        implicit none (type, external)

        ! Arguments
        integer, intent(in)        :: restart_header_status
        integer(int64), intent(out) :: nelmnt_a
        real(real64), intent(in)   :: closed_shell_breit_shift

        ! Local declarations
        real(real64)          :: densty
        ! Transfer parameters around and print out on node-0.
        !  .Temporary variables (ending with "_a") are used to store the
        !   sums over nodes;
        !  .Outputs to stream 24 are done on node-0;
        ! The following parameters are accumulated in setham which will not
        ! contain the correct values in restart mode. And thus is skipped.
        if (restart_header_status /= 0) then                    ! non-restart mode
            write (summary_output_unit, '("CUTOFF set to ",1p,d17.10)') block_matrix_element_cutoff
            write (summary_output_unit, '("Dirac-Coulomb one-e radial integrals:",1i8)') block_coulomb_one_body_integral_count
            write (summary_output_unit, '("One-e angular integrals that exceed CUTOFF: ",1i8)') &
                block_coulomb_one_body_coefficient_count
            write (summary_output_unit, '("Coulomb two-e radial integrals: ",1i8)') block_slater_integral_count
            write (summary_output_unit, '("Two-e angular integrals that exceed CUTOFF: ",1i8)') &
                block_coulomb_two_body_coefficient_count
            if (include_transverse_interaction) then
                write (summary_output_unit, '("Transverse two-e radial integrals: "/,6i8)') block_breit_integral_counts
                write (summary_output_unit, '("MCBP coefficients that exceed CUTOFF: ",1i8)') block_breit_coefficient_count
                write (summary_output_unit, '("Core coefficients that exceed CUTOFF: ",1i8)') block_breit_core_coefficient_count
            end if
            if (include_vacuum_polarization) write (summary_output_unit, '("Vacuum polarisation integrals: ",1i8)') &
                block_vacuum_polarization_integral_count
            if (include_normal_mass_shift) write (summary_output_unit, '("Kinetic energy integrals: ",1i8)') &
                block_kinetic_integral_count
            if (include_specific_mass_shift) write (summary_output_unit, '("Vinti integrals: ",1i8)') block_vinti_integral_count
        else
            write (summary_output_unit, *) "Restart mode --- no report on radial integrals"
        end if                                      !(irestart .NE. 0) ! non-restart mode

        ! ELSTO, EAV are not only for print-out, but also used later.
        ! density of the Hamiltonian matrix is only for print-out.
        ! At this place Hamiltonian matrix elements (EMT on each node) do
        !   _not_ contain ELSTO. ELSTO will be added to the total energy
        !  later with EAV.

        nelmnt_a = block_stored_element_count

        densty = dble(block_stored_element_count)/dble((block_assembly_csf_count*(block_assembly_csf_count + 1))/2)

        hamiltonian_energy_shift = hamiltonian_energy_shift/dble(block_assembly_csf_count) + closed_shell_breit_shift

        write (summary_output_unit, '("Elements that exceed CUTOFF in the lower"," triangle of the H matrix: ",1i8)') &
            block_stored_element_count
        write (summary_output_unit, '("Density of the H(amiltonian) matrix: ",1p,d22.15)') densty
        write (summary_output_unit, *)
        write (summary_output_unit, '("Average energy = ",1p,d19.12," Hartrees.")') hamiltonian_energy_shift
        write (summary_output_unit, *) hamiltonian_energy_shift, closed_shell_breit_shift, block_assembly_csf_count

        return
    end subroutine assemble_ci_matrix_from_restart
end module assemble_ci_matrix_from_restart_interface

module read_state_mixing_coefficients_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_state_mixing_coefficients
contains

    ! Open, check, load data from and close the rscf.mix file.
    !
    ! Written by Farid A. Parpia. Last revision: 25 Dec 1992
    subroutine read_state_mixing_coefficients(name, inpci, iblk)
        ! Module dependencies
        use calculation_state_files, only: mixing_input_path
        use calculation_file_units, only: mixing_input_unit
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc
        use atomic_calculation_state, only: electron_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use eigenstate_expansions, only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
        use selected_level_indices, only: active_level_indices, active_level_count
        use level_symmetry_metadata, only: level_parity_signs, level_two_j_plus_one
        use transition_basis_blocks, only: mixing_file_electron_count, mixing_file_csf_count, mixing_file_orbital_count, &
            mixing_file_level_count, mixing_file_coefficient_count, mixing_file_block_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: inpci
        integer, intent(in) :: iblk
        character(len=*), intent(in) :: name

        ! Local declarations
        integer :: j, io_status, nb, iatjp, iaspa, i
        character :: mixing_file_signature*6

        if (iblk == 1) then
            open(newunit=mixing_input_unit, file=mixing_input_path(name, inpci == 0), &
                form="UNFORMATTED", status="OLD", position="asis")

            read (mixing_input_unit, iostat=io_status) mixing_file_signature
            if (io_status/=0 .or. mixing_file_signature/="G92MIX") then
                write (output_unit, *)                                              &
                    "File", mixing_input_unit, "Not a GRASP MIXing Coefficients File;"
                close(mixing_input_unit)
                stop
            end if
            read (mixing_input_unit) mixing_file_electron_count, mixing_file_csf_count, mixing_file_orbital_count, &
                mixing_file_level_count, mixing_file_coefficient_count, mixing_file_block_count
            if (electron_count/=mixing_file_electron_count .or. orbital_count/=mixing_file_orbital_count) then
                ! :    (NCF .NE. NCFT) .OR.
                write (output_unit, *) "File", mixing_input_unit, "is not an"
                write (output_unit, *) " appropriate to Coefficients file"
                close(mixing_input_unit)
                stop
            end if
        end if
        ! Load data from the  rscf.mix  file
        write (output_unit, *) "Loading MIXing Coefficients File ..."
        read (mixing_input_unit) nb, csf_count, active_level_count, iatjp, iaspa
        call alloc (shifted_eigenvalues, active_level_count, "shifted_eigenvalues", "read_state_mixing_coefficients")
        call alloc (csf_mixing_coefficients, csf_count*active_level_count, "csf_mixing_coefficients", &
            "read_state_mixing_coefficients")
        call alloc (active_level_indices, active_level_count, "active_level_indices", "read_state_mixing_coefficients")
        call alloc (level_two_j_plus_one, active_level_count, "level_two_j_plus_one", "read_state_mixing_coefficients")
        call alloc (level_parity_signs, active_level_count, "level_parity_signs", "read_state_mixing_coefficients")
        ! These arrays are deallocated in mcp
        read (mixing_input_unit) (active_level_indices(i),i=1,active_level_count)
        level_two_j_plus_one(:active_level_count) = iatjp
        level_parity_signs(:active_level_count) = iaspa
        read (mixing_input_unit) hamiltonian_energy_shift, (shifted_eigenvalues(i),i=1,active_level_count)
        read (mixing_input_unit) ((csf_mixing_coefficients(i + (j - 1)*csf_count),i=1,csf_count),j=1,active_level_count)
        write (output_unit, *) " ... load complete;"
        ! Close the  rscf.mix  file
        return
    end subroutine read_state_mixing_coefficients
end module read_state_mixing_coefficients_interface

! Open the mixing output file and write its program-specific header.
module open_mixing_output_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mixing_output
contains

    ! Opens the  .mix  file on stream 25; writes a header to this file.
    !
    ! Written by Farid A. Parpia. Last revision: 24 Dec 1992
    ! Modified by Xinghong He. Last revision: 13 Jul 1998
    subroutine open_mixing_output(name)
        ! Module dependencies
        use calculation_file_units, only: mixing_output_unit
        use atomic_calculation_state, only: electron_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use hamiltonian_symmetry_blocks, only: symmetry_block_count
        use standard_io_units, only: stderr_unit
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        character, intent(in) :: name*(*)

        ! Local declarations
        integer :: open_status
        ! !...Wants nblock only
        call open_calculation_file (mixing_output_unit, name, "UNFORMATTED", "NEW", open_status)
        if (open_status /= 0) then
            write (stderr_unit, *) "Error when opening ", name(1:len_trim(name))
            stop
        end if
        ! Write the file header
        write (mixing_output_unit) "G92MIX"
        write (mixing_output_unit) electron_count, csf_count, orbital_count, 0, 0, symmetry_block_count
        ! ...The above record will be overidden in matrix.f
        !    with the final form of

        return
    end subroutine open_mixing_output
end module open_mixing_output_interface

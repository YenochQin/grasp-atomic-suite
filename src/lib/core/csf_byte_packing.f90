! Store a subshell value in a signed byte after checking its range and the destination index.
module csf_byte_packing
    implicit none (type, external)

    ! Module interface
    private
    public :: pack_csf_subshell_byte
contains
    subroutine pack_csf_subshell_byte(unpacked_value, byte_index, packed_values)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int8
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: unpacked_value
        integer, intent(in) :: byte_index
        integer(int8), intent(inout) :: packed_values(*)
        if (abs(unpacked_value) > 127) then
            write (stderr_unit, *) "pack_csf_subshell_byte: value out of range."
            stop
        end if
        if (byte_index <= 0) then
            write (stderr_unit, *) "pack_csf_subshell_byte: byte index ", byte_index, " less than 1"
            stop
        end if
        packed_values(byte_index) = unpacked_value
        return
    end subroutine pack_csf_subshell_byte
end module csf_byte_packing

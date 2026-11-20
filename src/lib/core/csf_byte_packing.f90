! Store a subshell value in a signed byte after checking its range and the destination index.
module csf_byte_packing_impl
    implicit none (type, external)

    ! Module interface
    private
    public :: pack_impl
contains
    subroutine pack_impl(iunpkd, isubsh, ipackd)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int8
        use standard_io_units, only: stderr_unit
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: iunpkd
        integer, intent(in) :: isubsh
        integer(int8), intent(inout) :: ipackd(*)
        if (abs(iunpkd) > 127) then
            write (stderr_unit, *) "PACK: Argument IUNPKD out of range."
            stop
        end if
        if (isubsh <= 0) then
            write (stderr_unit, *) "PACK: ISUBSH=", isubsh, " less than 1"
            stop
        end if
        ipackd(isubsh) = iunpkd
        return
    end subroutine pack_impl
end module csf_byte_packing_impl

! Export PACK as a module procedure to keep it distinct from the Fortran intrinsic in callers.
module csf_byte_packing
    ! Module dependencies
    use csf_byte_packing_impl, only: pack => pack_impl
    implicit none (type, external)

    ! Module interface
    private
    public :: pack
end module csf_byte_packing

module right_csf_two_j_plus_one_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_two_j_plus_one
contains

    ! ITJPOR is the value of 2J+1 for CSF number ICSF.
    !
    ! Written by Farid A. Parpia. Last revision: 23 Dec 1992
    ! Modified by G. Gaigalas. May 2011
    integer function right_csf_two_j_plus_one (csf_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use csf_coupling_data,        only: left_angular_codes, csf_total_j_slot
        use orbital_csf_basis,         only: left_csf_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_index

        if (csf_index>=1 .and. csf_index<=left_csf_count) then
            right_csf_two_j_plus_one = left_angular_codes(csf_total_j_slot,csf_index)
        else
            write (output_unit, *) "right_csf_two_j_plus_one: Argument csf_index is out of range."
            stop
        end if
        return
    end function right_csf_two_j_plus_one
end module right_csf_two_j_plus_one_interface

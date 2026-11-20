module right_csf_two_j_plus_one_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_two_j_plus_one
contains

    !***********************************************************************
    !                                                                      *
    integer function right_csf_two_j_plus_one (csf_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !   Modified by G. Gaigalas                                 May 2011   *
        !                                                                      *
        !***********************************************************************
        use grasp_parameters, only: max_orbitals
        use csf_coupling_data,        only: left_coupling_codes
        use orbital_csf_basis,         only: left_csf_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_index
        !-----------------------------------------------
        !
        if (csf_index>=1 .and. csf_index<=left_csf_count) then
            right_csf_two_j_plus_one = left_coupling_codes(max_orbitals,csf_index)
            if (right_csf_two_j_plus_one > 127) right_csf_two_j_plus_one = 256 - right_csf_two_j_plus_one
            right_csf_two_j_plus_one = iabs (right_csf_two_j_plus_one)
        else
            write (output_unit, *) "ITJPOR: Argument ICSF is out of range."
            stop
        end if
        !
        return
    end function right_csf_two_j_plus_one
end module right_csf_two_j_plus_one_interface

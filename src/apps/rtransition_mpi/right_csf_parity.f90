module right_csf_parity_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_parity
contains

    !***********************************************************************
    !                                                                      *
    integer function right_csf_parity (csf_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   ISPARR is the value of P for CSF number ICSF.                      *
        !                                                                      *
        !   Call(s) to: [LIB92]: IUNPCK.                                       *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !   Modified by G. Gaigalas                                 May 2011   *
        !                                                                      *
        !***********************************************************************
        use csf_coupling_data,        only: left_angular_codes, csf_parity_slot
        use orbital_csf_basis,         only: left_csf_count
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: csf_index
        !
        !
        if (csf_index>=1 .and. csf_index<=left_csf_count) then
            right_csf_parity = left_angular_codes(csf_parity_slot,csf_index)

        else
            write (output_unit, *) "right_csf_parity: Argument csf_index is out of range."
            stop
        end if
        !
        return
    end function right_csf_parity
end module right_csf_parity_interface

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
        !
        !
        if (csf_index>=1 .and. csf_index<=left_csf_count) then
            right_csf_parity = left_coupling_codes(max_orbitals,csf_index)

            if (right_csf_parity > 127) right_csf_parity = right_csf_parity - 256
            right_csf_parity = sign(1,right_csf_parity)
        else
            write (output_unit, *) "ISPARR: Argument ICSF is out of range."
            stop
        end if
        !
        return
    end function right_csf_parity
end module right_csf_parity_interface

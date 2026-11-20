module right_csf_intermediate_coupling_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_intermediate_coupling
contains

    !***********************************************************************
    !                                                                      *
    integer function right_csf_intermediate_coupling (loc, csf_index)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   JCUPR is the 2J+1 value of the LOCth nontrivial intermediate ang-  *
        !   ular momentum in CSF  ICSF.                                        *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !   Modified by G. Gaigalas                                 May 2011   *
        !                                                                      *
        !***********************************************************************
        use orbital_csf_basis,         only: left_orbital_count, left_csf_count
        use csf_coupling_data,        only: left_coupling_codes
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: loc
        integer, intent(in) :: csf_index
        !-----------------------------------------------
        !
        if (loc>=1 .and. loc<=left_orbital_count-1) then
            if (csf_index>=1 .and. csf_index<=left_csf_count) then
                right_csf_intermediate_coupling = left_coupling_codes(loc,csf_index)
            else
                write (output_unit, *) "JCUPR: Argument ICSF is out of range."
                stop
            end if
        else
            write (output_unit, *) "JCUPR: Argument LOC is out of range."
            stop
        end if
        !
        return
    end function right_csf_intermediate_coupling
end module right_csf_intermediate_coupling_interface

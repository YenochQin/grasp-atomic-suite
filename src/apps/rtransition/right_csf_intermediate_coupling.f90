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
    use grasp_parameters, only: max_orbitals
    use orbital_csf_basis,         only: nwr, ncfr
    use csf_coupling_data,        only: jcupar
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: loc
    integer, intent(in) :: csf_index
    !-----------------------------------------------
    !
    if (loc>=1 .and. loc<=nwr-1) then
        if (csf_index>=1 .and. csf_index<=ncfr) then
            right_csf_intermediate_coupling = jcupar(loc,csf_index)
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

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
    use grasp_parameters, only: max_orbitals
    use csf_coupling_data,        only: jcupar
    use standard_io_units,      only: stderr_unit
    use orbital_csf_basis,         only: ncfr
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: csf_index
    !
    !
    if (csf_index>=1 .and. csf_index<=ncfr) then
        right_csf_parity = jcupar(max_orbitals,csf_index)

        if (right_csf_parity > 127) right_csf_parity = right_csf_parity - 256
        right_csf_parity = sign(1,right_csf_parity)
    else
        write (output_unit, *) "ISPARR: Argument ICSF is out of range."
        stop
    end if
    !
    return
end function right_csf_parity

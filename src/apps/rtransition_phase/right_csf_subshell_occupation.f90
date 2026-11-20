module right_csf_subshell_occupation_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_subshell_occupation
contains

    !***********************************************************************
    !                                                                      *
    integer function right_csf_subshell_occupation (subshell_index, csf_index)
        ! Module dependencies
        !                                                                      *
        !   IQR is the occupation of subshell ISUBSH in CSF  ICSF.             *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !   Modified by G. Gaigalas                                 May 2011   *
        !                                                                      *
        !***********************************************************************
        use orbital_csf_basis,         only: csf_subshell_occupations
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: subshell_index
        integer, intent(in) :: csf_index
        !-----------------------------------------------
        !
        !
        right_csf_subshell_occupation = csf_subshell_occupations(subshell_index,csf_index)
        !
        return
    end function right_csf_subshell_occupation
end module right_csf_subshell_occupation_interface

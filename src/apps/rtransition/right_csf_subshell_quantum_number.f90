module right_csf_subshell_quantum_number_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: right_csf_subshell_quantum_number
contains

    !***********************************************************************
    !                                                                      *
    integer function right_csf_subshell_quantum_number (iwhich, subshell_index, csf_index)
        ! Module dependencies
        !                                                                      *
        !   JQSR is a subshell quantum number for subshell ISUBSH in configu-  *
        !   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
        !   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
        !   Modified by G. Gaigalas                                 May 2011   *
        !                                                                      *
        !***********************************************************************
        use csf_coupling_data,        only: left_subshell_quantum_codes
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: iwhich
        integer, intent(in) :: subshell_index
        integer, intent(in) :: csf_index
        !-----------------------------------------------
        !
        right_csf_subshell_quantum_number = left_subshell_quantum_codes(subshell_index,iwhich,csf_index)
        !
        return
    end function right_csf_subshell_quantum_number
end module right_csf_subshell_quantum_number_interface

module right_csf_subshell_quantum_number_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function right_csf_subshell_quantum_number (iwhich, subshell_index, csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: iwhich
            integer, intent(in) :: subshell_index
            integer, intent(in) :: csf_index
        end function right_csf_subshell_quantum_number
    end interface
end module right_csf_subshell_quantum_number_interface

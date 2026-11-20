module copy_state_csf_data_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine copy_state_csf_data (n)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n
        end subroutine copy_state_csf_data
    end interface
end module copy_state_csf_data_interface

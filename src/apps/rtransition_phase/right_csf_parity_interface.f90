module right_csf_parity_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function right_csf_parity (csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_index
        end function right_csf_parity
    end interface
end module right_csf_parity_interface

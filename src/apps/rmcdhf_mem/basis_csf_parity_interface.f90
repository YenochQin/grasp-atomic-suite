module basis_csf_parity_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        integer function basis_csf_parity (csf_index)
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_index
        end function basis_csf_parity
    end interface
end module basis_csf_parity_interface

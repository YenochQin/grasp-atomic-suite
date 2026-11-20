module write_csf_weights_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_csf_weights (write_standard_output)
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: write_standard_output
        end subroutine write_csf_weights
    end interface
end module write_csf_weights_interface

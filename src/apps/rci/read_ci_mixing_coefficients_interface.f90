module read_ci_mixing_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_ci_mixing_coefficients (block_labels)
            implicit none (type, external)

            ! Arguments
            character(len=8), intent(in):: block_labels(*)
        end subroutine read_ci_mixing_coefficients
    end interface
end module read_ci_mixing_coefficients_interface

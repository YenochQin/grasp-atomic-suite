module read_al_level_weights_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_al_level_weights (csf_count, level_weights)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_count
            real(real64), intent(inout) :: level_weights(csf_count)
        end subroutine read_al_level_weights
    end interface
end module read_al_level_weights_interface

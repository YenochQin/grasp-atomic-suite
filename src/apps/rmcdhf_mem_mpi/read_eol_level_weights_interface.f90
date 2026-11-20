module read_eol_level_weights_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_eol_level_weights (ndef, ncmin, level_weights)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer :: ndef
            !VAST...Dummy argument NDEF is not referenced in this routine.
            integer, intent(in) :: ncmin
            real(real64), intent(inout) :: level_weights(ncmin)
        end subroutine read_eol_level_weights
    end interface
end module read_eol_level_weights_interface

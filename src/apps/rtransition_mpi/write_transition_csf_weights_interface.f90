module write_transition_csf_weights_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_transition_csf_weights (asfa,asfb,lev1,lev2)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: asfa, asfb
            integer, intent(in) :: lev1, lev2
        end subroutine write_transition_csf_weights
    end interface
end module write_transition_csf_weights_interface

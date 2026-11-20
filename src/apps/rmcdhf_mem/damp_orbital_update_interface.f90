module damp_orbital_update_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine damp_orbital_update (orbital_index, orbital_sign_reversal_count, orbital_damping_factor)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: orbital_index
            integer, intent(inout) :: orbital_sign_reversal_count
            real(real64), intent(in) :: orbital_damping_factor
        end subroutine damp_orbital_update
    end interface
end module damp_orbital_update_interface

module second_spin_angular_factor_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine second_spin_angular_factor (ipar,i1,i2,apart)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: ipar, i1, i2
            real(real64), intent(out) :: apart
        end subroutine second_spin_angular_factor
    end interface
end module second_spin_angular_factor_interface

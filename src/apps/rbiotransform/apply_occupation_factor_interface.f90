module apply_occupation_factor_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine apply_occupation_factor (input_mixing_coefficients, csf_count, level_count, i, l, const, &
            & transformed_mixing_coefficients, ntestg)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_count
            integer, intent(in) :: level_count
            real(real64), intent(in) :: input_mixing_coefficients(csf_count, level_count)
            integer, intent(in) :: i
            integer, intent(in) :: l
            real(real64), intent(in) :: const
            real(real64), intent(inout) :: transformed_mixing_coefficients(csf_count, level_count)
            integer, intent(in) :: ntestg
        end subroutine apply_occupation_factor
    end interface
end module apply_occupation_factor_interface

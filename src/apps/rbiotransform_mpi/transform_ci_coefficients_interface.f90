module transform_ci_coefficients_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine transform_ci_coefficients (input_mixing_coefficients,csf_count,level_count,l,subshell_count,t,nin,ntestg, &
            & transformed_mixing_coefficients,scr)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: csf_count
            integer, intent(in) :: level_count
            real(real64) :: input_mixing_coefficients(csf_count, level_count)
            integer, intent(in) :: l
            integer, intent(in) :: subshell_count
            real(real64), intent(in) :: t(subshell_count, subshell_count)
            integer, intent(in) :: nin
            integer, intent(in) :: ntestg
            real(real64) :: transformed_mixing_coefficients(csf_count, level_count)
            real(real64) :: scr(csf_count, level_count)
        end subroutine transform_ci_coefficients
    end interface
end module transform_ci_coefficients_interface

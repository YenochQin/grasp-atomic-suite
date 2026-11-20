module write_transition_probability_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine write_transition_probability (asfa, asfb, i, j, omega, factor, lines, lsame)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            logical, intent(in) :: lsame
            real(real64), intent(in) :: asfa
            real(real64), intent(in) :: asfb
            integer :: i
            integer, intent(in) :: j
            real(real64), intent(in) :: omega
            real(real64), intent(in) :: factor
            integer, intent(inout) :: lines
        end subroutine write_transition_probability
    end interface
end module write_transition_probability_interface

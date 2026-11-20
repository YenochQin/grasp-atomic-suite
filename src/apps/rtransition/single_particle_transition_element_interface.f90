module single_particle_transition_element_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine single_particle_transition_element (first_orbital, second_orbital, hcoul, hbab, hmag)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: first_orbital
            integer, intent(in) :: second_orbital
            real(real64), intent(out) :: hcoul
            real(real64), intent(out) :: hbab
            real(real64), intent(out) :: hmag
        end subroutine single_particle_transition_element
    end interface
end module single_particle_transition_element_interface

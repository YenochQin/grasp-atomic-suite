module check_orbital_damping_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine check_orbital_damping (ipr, orbital_index, ed1, ed2)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(inout) :: ipr
            integer, intent(in) :: orbital_index
            real(real64), intent(inout) :: ed1
            real(real64), intent(inout) :: ed2
        end subroutine check_orbital_damping
    end interface
end module check_orbital_damping_interface

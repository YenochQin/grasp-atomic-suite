module multiply_restart_hamiltonian_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine multiply_restart_hamiltonian (n, m, b, c)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: n
            integer, intent(in) :: m
            real(real64), intent(in) :: b(n, m)
            real(real64), intent(inout) :: c(n, m)
        end subroutine multiply_restart_hamiltonian
    end interface
end module multiply_restart_hamiltonian_interface

module update_scf_mixing_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine update_scf_mixing (sum)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: sum
        end subroutine update_scf_mixing
    end interface
end module update_scf_mixing_interface

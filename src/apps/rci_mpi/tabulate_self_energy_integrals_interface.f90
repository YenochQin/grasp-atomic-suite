module tabulate_self_energy_integrals_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine tabulate_self_energy_integrals (slfint)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_orbitals
            implicit none (type, external)

            ! Arguments
            real(real64), intent(out) :: slfint(max_orbitals)
        end subroutine tabulate_self_energy_integrals
    end interface
end module tabulate_self_energy_integrals_interface

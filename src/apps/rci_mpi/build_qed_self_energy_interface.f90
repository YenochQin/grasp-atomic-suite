module build_qed_self_energy_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine build_qed_self_energy (jstate, slfint, ucf)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: jstate
            real(real64), intent(out) :: slfint(max_orbitals)
            real(real64), intent(out) :: ucf(:)
        end subroutine build_qed_self_energy
    end interface
end module build_qed_self_energy_interface

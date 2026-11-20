module read_specific_mass_shift_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine read_specific_mass_shift_mpi (vint,vint2)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
        end subroutine read_specific_mass_shift_mpi
    end interface
end module read_specific_mass_shift_mpi_interface

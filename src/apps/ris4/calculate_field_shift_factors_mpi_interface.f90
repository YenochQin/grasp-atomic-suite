module calculate_field_shift_factors_mpi_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        subroutine calculate_field_shift_factors_mpi (doit,dint1,dint2,dint3,dint4,dint5,dint6,dint7)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters, only: max_orbitals
            implicit none (type, external)

            ! Arguments
            integer, intent(in) :: doit
            real(real64), intent(in) :: dint1(max_orbitals, max_orbitals), dint2(max_orbitals, max_orbitals), &
                dint3(max_orbitals, max_orbitals), dint4(max_orbitals, max_orbitals), &
                dint5(max_orbitals, max_orbitals), dint6(max_orbitals, max_orbitals), &
                dint7(max_orbitals, max_orbitals)
        end subroutine calculate_field_shift_factors_mpi
    end interface
end module calculate_field_shift_factors_mpi_interface

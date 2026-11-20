module calculate_seltzer_field_shift_factors_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine calculate_seltzer_field_shift_factors (doit,dint1,dint2,dint3,dint4,dint5,dint6,    &
            dint7,dint1vec,dens1vec,nrnuc)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_orbitals, max_radial_points
            use selected_level_indices,          only : nvec
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: dint1(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint2(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint3(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint4(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint5(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint6(max_orbitals, max_orbitals)
            real(real64), intent(in) :: dint7(max_orbitals, max_orbitals)
            integer, intent(in) :: nrnuc
            real(real64), intent(out) :: dens1vec(nvec, nrnuc)
            real(real64), intent(in) :: dint1vec(max_orbitals, max_orbitals, nrnuc)
            integer, intent(in) :: doit
        end subroutine calculate_seltzer_field_shift_factors
    end interface
end module calculate_seltzer_field_shift_factors_interface

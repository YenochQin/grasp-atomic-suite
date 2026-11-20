module calculate_specific_mass_shift_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine calculate_specific_mass_shift (doit,vint,vint2)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_orbitals
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
            integer, intent(in) :: doit
        end subroutine calculate_specific_mass_shift
    end interface
end module calculate_specific_mass_shift_interface

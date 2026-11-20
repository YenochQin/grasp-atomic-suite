module read_specific_mass_shift_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine read_specific_mass_shift (vint,vint2)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,   only: max_orbitals
            implicit none (type, external)

            ! Arguments
            real(real64), intent(in) :: vint(max_orbitals, max_orbitals), vint2(max_orbitals, max_orbitals)
        end subroutine read_specific_mass_shift
    end interface
end module read_specific_mass_shift_interface

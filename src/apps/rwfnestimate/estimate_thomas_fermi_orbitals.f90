! Solve for initial radial orbitals in the Thomas–Fermi potential.
module estimate_thomas_fermi_orbitals_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: estimate_thomas_fermi_orbitals
contains

    !   This  subroutine  is  used  to  produce  estimates  of  the wave
    !   functions by use of the Thomas-Fermi approximation to the direct
    !   potential. SUBROUTINE  SOLVH  is used  to obtain the radial wave
    !   functions.
    !
    ! Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992
    subroutine estimate_thomas_fermi_orbitals(matched_orbital_indices, matched_subshell_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use initial_orbital_selection, only: orbital_initialized
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: max_orbitals, orbital_origin_parameters
        use initial_orbital_source, only: orbital_source_description

        ! Procedure interfaces
        use solve_homogeneous_radial_equation_interface, only: solve_homogeneous_radial_equation
        implicit none (type, external)

        ! Arguments
        !
        !   G l o b a l   P a r a m e t e r s
        integer , intent(in) :: matched_subshell_count
        integer , intent(in) :: matched_orbital_indices(max_orbitals)

        ! Local declarations
        integer :: matched_subshell_index, orbital_index
        logical :: solve_failed

        do matched_subshell_index = 1, matched_subshell_count
            orbital_index = matched_orbital_indices(matched_subshell_index)
            if (orbital_initialized(orbital_index)) cycle
            ! Estimate the leading coefficient of the series expansion for
            ! the orbital
            orbital_origin_parameters(orbital_index) = 10.0_real64
            ! Calculate radial wavefunctions
            call solve_homogeneous_radial_equation (orbital_index, solve_failed)
            ! Message if SOLVH did not converge; reset SET(LOC) otherwise
            if (solve_failed) then
                write (*, '(/,"TFWAVE: Unable to compute radial"/,/," wavefunction for ",i2,a2," subshell;")') &
                    orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
            else
                orbital_initialized(orbital_index) = .true.
                ! SOURCE(LOC) = 'Thomas-Fermi estimate'
                orbital_source_description(orbital_index) = "T-F"
            end if
        end do
        return
        return
    end subroutine estimate_thomas_fermi_orbitals
end module estimate_thomas_fermi_orbitals_interface

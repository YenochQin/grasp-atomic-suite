module write_transition_level_energies_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_transition_level_energies
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_transition_level_energies(hamiltonian_energy_shift, e, jtot, ipar, level_indices, level_count, mode, k)
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        !                                                                      *
        !   This  subroutine prints  energy levels, splittings, and energies   *
        !   relative to the lowest in  Hartrees, Kaysers, and  eV, using the   *
        !   reduced mass corrected value for the Rydberg. If MODE is 0, only   *
        !   the eigenenergies are printed. If  MODE  is 1, the eigenenergies   *
        !   and separations are printed. If MODE is 2, the eigenenergies and   *
        !   energies relative to level 1 are printed. If MODE is 3, the eig-   *
        !   enenergies,  separations,  and energies relative to level  1 are   *
        !   printed.                                                           *
        !                                          Last updated: 15 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: hartree_to_wavenumber_cm, hartree_to_ev
        use angular_symmetry_labels, only: format_j_label, parity_labels
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: level_count
        integer  :: mode
        integer , intent(in) :: k
        real(real64) , intent(in) :: hamiltonian_energy_shift
        integer , intent(in) :: jtot(level_count)
        integer , intent(in) :: ipar(level_count)
        integer , intent(in) :: level_indices(level_count)
        real(real64) , intent(in) :: e(level_count)

        ! Local declarations
        integer :: j, i, ip
        real(real64) :: eau, ecm, eev
        !-----------------------------------------------
        !
        !   Always print the eigenenergies
        !
        if (k == 1) then
            write (summary_output_unit, '("Eigenenergies for the initial state list")')
        end if
        if (k == 2) then
            write (summary_output_unit, '("Eigenenergies for the final state list")')
        end if
        write (summary_output_unit, '("Level  J Parity",10x,"Hartrees",18x,"Kaysers")')
        do j = 1, level_count
            i = level_indices(j)
            eau = e(j) + hamiltonian_energy_shift
            ecm = eau*hartree_to_wavenumber_cm
            eev = eau*hartree_to_ev
            ip = (ipar(j)+3)/2
            write (summary_output_unit, '(1i3,4x,2a,1p,2d25.15)') i, format_j_label(jtot(j)), parity_labels(ip), eau, ecm
        end do
        !
        return
        !
        return
        !
    end subroutine write_transition_level_energies
end module write_transition_level_energies_interface

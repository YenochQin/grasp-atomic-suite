module write_level_energies_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_level_energies
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_level_energies(hamiltonian_energy_shift, e, jtot, ipar, level_indices, level_count, mode)
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
        use atomic_calculation_state, only: hartree_to_ev, hartree_to_wavenumber_cm
        use angular_symmetry_labels, only: format_j_label, parity_labels
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: jtot
        integer, intent(in) :: ipar
        integer, intent(in) :: level_count
        integer, intent(in) :: mode
        real(real64), intent(in) :: hamiltonian_energy_shift
        integer, intent(in) :: level_indices(level_count)
        real(real64), intent(in) :: e(level_count)

        ! Local declarations
        integer      :: j, i, ip
        real(real64) :: eau, ecm, eev
        !-----------------------------------------------
        !
        !   Always print the eigenenergies
        !
        write (summary_output_unit, '(/,"Eigenenergies:")')
        write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
        do j = 1, level_count
            i = level_indices(j)
            eau = e(j) + hamiltonian_energy_shift
            ecm = eau*hartree_to_wavenumber_cm
            eev = eau*hartree_to_ev
            ip = (ipar + 3)/2
            write (summary_output_unit, '(1i3,2x,2a,1p,3d22.14)') i, format_j_label(jtot), parity_labels(ip), eau, &
                & ecm, eev
        end do
        !
        if (level_count > 1) then
            !
            !   Energy separations
            !
            if (mode==1 .or. mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do j = 2, level_count
                    i = level_indices(j)
                    eau = e(j) - e(j-1)
                    ecm = eau*hartree_to_wavenumber_cm
                    eev = eau*hartree_to_ev
                    ip = (ipar + 3)/2
                    write (summary_output_unit, '(1i3,2x,2a,1p,3d22.14)') i, format_j_label(jtot), parity_labels(ip), &
                        & eau, ecm, eev
                end do
            end if
            !
            !   Energies relative to level 1
            !
            if (mode==2 .or. mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to lowest level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do j = 2, level_count
                    i = level_indices(j)
                    eau = e(j) - e(1)
                    ecm = eau*hartree_to_wavenumber_cm
                    eev = eau*hartree_to_ev
                    ip = (ipar + 3)/2
                    write (summary_output_unit, '(1i3,2x,2a,1p,3d22.14)') i, format_j_label(jtot), parity_labels(ip), &
                        & eau, ecm, eev
                end do
            end if
            !
        end if
        !
        return
        !
        return
        !
    end subroutine write_level_energies
end module write_level_energies_interface

module write_hyperfine_energies_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_hyperfine_energies
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_hyperfine_energies(hamiltonian_energy_shift, shifted_eigenvalues, level_two_j_plus_one, level_parity_signs, &
                & level_indices, level_count, energy_separation_mode)
        ! Module dependencies
        use calculation_file_units, only: hyperfine_output_unit, summary_output_unit
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
        use angular_symmetry_labels, only: right_justified_j_labels, parity_labels
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: level_count
        integer, intent(in) :: energy_separation_mode
        real(real64), intent(in) :: hamiltonian_energy_shift
        integer, intent(in) :: level_two_j_plus_one(level_count)
        integer, intent(in) :: level_parity_signs(level_count)
        integer, intent(in) :: level_indices(level_count)
        real(real64), intent(in) :: shifted_eigenvalues(level_count)

        ! Local declarations
        integer      :: level_index, level_label, parity_label_index
        real(real64) :: energy_hartree, energy_wavenumber, energy_ev
        !-----------------------------------------------
        !
        !   Always print the eigenenergies
        !
        write (summary_output_unit, '(/,"Eigenenergies:")')
        write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
        write (hyperfine_output_unit, '(/,"Eigenenergies:")')
        write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
        do level_index = 1, level_count
            level_label = level_indices(level_index)
            energy_hartree = shifted_eigenvalues(level_index) + hamiltonian_energy_shift
            energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
            energy_ev = energy_hartree*hartree_to_ev
            parity_label_index = (level_parity_signs(level_index)+3)/2
            write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), parity_labels(parity_label_index), energy_hartree, &
                & energy_wavenumber, energy_ev
            write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), parity_labels(parity_label_index), energy_hartree, &
                & energy_wavenumber, energy_ev
        end do
        !
        if (level_count > 1) then
            !
            !   Energy separations
            !
            if (energy_separation_mode==1 .or. energy_separation_mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                write (hyperfine_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
                write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do level_index = 2, level_count
                    level_label = level_indices(level_index)
                    energy_hartree = shifted_eigenvalues(level_index) - shifted_eigenvalues(level_index-1)
                    energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
                    energy_ev = energy_hartree*hartree_to_ev
                    parity_label_index = (level_parity_signs(level_index)+3)/2
                    write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), &
                        & parity_labels(parity_label_index), &
                        & energy_hartree, energy_wavenumber, energy_ev
                    write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), &
                        & parity_labels(parity_label_index), energy_hartree, energy_wavenumber, energy_ev
                end do
            end if
            !
            !   Energies relative to level 1
            !
            if (energy_separation_mode==2 .or. energy_separation_mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to lowest level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                write (hyperfine_output_unit, '(/,"Energy of each level relative to lowest level:")')
                write (hyperfine_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do level_index = 2, level_count
                    level_label = level_indices(level_index)
                    energy_hartree = shifted_eigenvalues(level_index) - shifted_eigenvalues(1)
                    energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
                    energy_ev = energy_hartree*hartree_to_ev
                    parity_label_index = (level_parity_signs(level_index)+3)/2
                    write (summary_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), &
                        & parity_labels(parity_label_index), &
                        & energy_hartree, energy_wavenumber, energy_ev
                    write (hyperfine_output_unit, '(1i3,2x,2a4,1p,3d22.14)') level_label, &
                & right_justified_j_labels(level_two_j_plus_one(level_index)), &
                        & parity_labels(parity_label_index), energy_hartree, energy_wavenumber, energy_ev
                end do
            end if
            !
        end if
        !
        return
        !
    end subroutine write_hyperfine_energies
end module write_hyperfine_energies_interface

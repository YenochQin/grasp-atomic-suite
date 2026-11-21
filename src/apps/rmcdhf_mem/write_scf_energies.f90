module write_scf_energies_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_scf_energies
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_scf_energies(shifted_level_energies,level_indices,level_count,energy_difference_mode)
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
        !   Midified by G. Gaigalas                              05 Feb 2017   *
        !      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
        !                                  JCUPA(NNNW*NCF)                     *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: hartree_to_ev, hartree_to_wavenumber_cm
        use angular_symmetry_labels, only: format_j_label, parity_labels
        use transition_block_indices, only: level_block_indices
        use block_mean_energies, only: mean_energy_by_block
        use level_symmetry_metadata,          only: packed_block_symmetries
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: level_count
        integer, intent(in) :: energy_difference_mode
        integer, intent(in) :: level_indices(level_count)
        real(real64), intent(in) :: shifted_level_energies(level_count)

        ! Local declarations
        integer :: selected_level_index, symmetry_block_index, level_number, parity_label_index, two_j_plus_one
        real(real64) :: hamiltonian_energy_shift, energy_hartree, energy_wavenumber, energy_ev
        !-----------------------------------------------
        !
        !
        !   Always print the eigenenergies
        !
        write (summary_output_unit, '(/,"Eigenenergies:")')
        write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
        do selected_level_index = 1, level_count
            symmetry_block_index = level_block_indices(selected_level_index)
            hamiltonian_energy_shift = mean_energy_by_block(symmetry_block_index)
            level_number = level_indices(selected_level_index)
            energy_hartree = shifted_level_energies(selected_level_index) + hamiltonian_energy_shift
            energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
            energy_ev = energy_hartree*hartree_to_ev
            !GG         IP = (IPAR(J)+3)/2
            two_j_plus_one = iabs(packed_block_symmetries(symmetry_block_index))
            if(packed_block_symmetries(symmetry_block_index) >= 0) then
                parity_label_index = 2
            else
                parity_label_index = 1
            end if
            write (summary_output_unit, '(1i3,4x,2a,1p,3d21.12)') level_number, format_j_label(two_j_plus_one), &
                & parity_labels(parity_label_index), energy_hartree, &
                & energy_wavenumber, energy_ev
            !GG         WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
        end do
        !
        if (level_count > 1) then
            !
            !   Energy separations
            !
            if (energy_difference_mode==1 .or. energy_difference_mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to immediately lower"," level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do selected_level_index = 2, level_count
                    level_number = level_indices(selected_level_index)
                    energy_hartree = shifted_level_energies(selected_level_index) - shifted_level_energies(selected_level_index-1)
                    energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
                    energy_ev = energy_hartree*hartree_to_ev
                    !GG               IP = (IPAR(J)+3)/2
                    symmetry_block_index = level_block_indices(selected_level_index)
                    two_j_plus_one = iabs(packed_block_symmetries(symmetry_block_index))
                    if(packed_block_symmetries(symmetry_block_index) >= 0) then
                        parity_label_index = 2
                    else
                        parity_label_index = 1
                    end if
                    write (summary_output_unit, '(1i3,4x,2a,1p,3d21.12)') level_number, &
                        & format_j_label(two_j_plus_one), parity_labels(parity_label_index), &
                        & energy_hartree, energy_wavenumber, energy_ev
                    !GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
                end do
            end if
            !
            !   Energies relative to level 1
            !
            if (energy_difference_mode==2 .or. energy_difference_mode==3) then
                write (summary_output_unit, '(/,"Energy of each level relative to lowest level:")')
                write (summary_output_unit, '(/,"Level  J Parity",7x,"Hartrees",14x,"Kaysers",16x,"eV"/)')
                do selected_level_index = 2, level_count
                    level_number = level_indices(selected_level_index)
                    energy_hartree = shifted_level_energies(selected_level_index) - shifted_level_energies(1)
                    energy_wavenumber = energy_hartree*hartree_to_wavenumber_cm
                    energy_ev = energy_hartree*hartree_to_ev
                    !GG               IP = (IPAR(J)+3)/2
                    symmetry_block_index = level_block_indices(selected_level_index)
                    two_j_plus_one = iabs(packed_block_symmetries(symmetry_block_index))
                    if(packed_block_symmetries(symmetry_block_index) >= 0) then
                        parity_label_index = 2
                    else
                        parity_label_index = 1
                    end if
                    write (summary_output_unit, '(1i3,4x,2a,1p,3d21.12)') level_number, &
                        & format_j_label(two_j_plus_one), parity_labels(parity_label_index), &
                        & energy_hartree, energy_wavenumber, energy_ev
                    !GG               WRITE (24, 302) I, LABJ(JTOT(J)), LABP(IP), EAU, ECM, EEV
                end do
            end if
            !
        end if
        !
        return
        !
        return
        !
    end subroutine write_scf_energies
end module write_scf_energies_interface

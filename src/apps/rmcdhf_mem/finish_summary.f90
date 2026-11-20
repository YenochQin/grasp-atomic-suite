module finish_summary_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: finish_summary
contains

    !***********************************************************************
    !                                                                      *
    subroutine finish_summary
        ! Module dependencies
        use calculation_file_units, only: summary_output_unit
        !                                                                      *
        !   Generates the last part of  rscf92.sum  (on stream 24).            *
        !                                                                      *
        !   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
        !   Midified by G. Gaigalas                              05 Feb 2017   *
        !      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
        !                                  JCUPA(NNNW*NCF)                     *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use atomic_calculation_state, only: requested_level_indices, selected_level_count
        use eigenstate_expansions, only: shifted_eigenvalues
        use orbital_csf_basis, only: orbital_count, orbital_energy_parameters, orbital_origin_exponents, orbital_principal_n, &
            orbital_symmetry_labels
        use scf_iteration_state, only: orbital_occupation_weights, orbital_self_consistency_errors
        !gg      USE syma_C
        use radial_wavefunctions, only: orbital_last_grid_index, orbital_origin_parameters, radial_large_components, &
            radial_small_components
        use weighted_radial_overlap, only: integrate_weighted_radial_overlap

        ! Procedure interfaces
        use write_scf_energies_interface, only: write_scf_energies
        use write_csf_weights_interface, only: write_csf_weights
        implicit none (type, external)

        ! Local declarations
        integer :: orbital_index, mode
        real(real64) :: wa, wb, wc, wd, we
        !-----------------------------------------------
        !
        !   Write out the orbital properties
        !
        write (summary_output_unit, &
            '(/,"Radial wavefunction summary:"/,/,67x,"Self"/,"Subshell",6x,"e",13x,"p0",5x,"gamma",5x,"P(2)",&
            &7x,"Q(2)",3x,"Consistency"," MTP"/)')
        do orbital_index = 1, orbital_count
            write (summary_output_unit, '(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,3(d11.3),i5)') &
                orbital_principal_n(orbital_index),orbital_symmetry_labels(orbital_index),orbital_energy_parameters( &
                & orbital_index), &
                & orbital_origin_parameters(orbital_index),orbital_origin_exponents(orbital_index),radial_large_components(2, &
                & orbital_index),radial_small_components(2,orbital_index), &
                orbital_self_consistency_errors(orbital_index), orbital_last_grid_index(orbital_index)
        end do
        !
        write (summary_output_unit, '(/18x,"-3",14x,"-1",29x,"2",14x,"4",5x,"Generalised"/"Subshell",4x,"<  r  >",8x,"<  r  >",8x,&
            &"<  r  >",8x,"<  r  >",8x,"<  r  >",6x,"occupation"/)')
        do orbital_index = 1, orbital_count
            wa = integrate_weighted_radial_overlap(orbital_index,orbital_index,-1)
            wb = integrate_weighted_radial_overlap(orbital_index,orbital_index,1)
            wc = integrate_weighted_radial_overlap(orbital_index,orbital_index,2)
            wd = integrate_weighted_radial_overlap (orbital_index,orbital_index, 4)
            we = 0.0_real64
            if (orbital_symmetry_labels(orbital_index) /= "s " .and. orbital_symmetry_labels(orbital_index) /= "p-") then
                we = integrate_weighted_radial_overlap(orbital_index,orbital_index,-3)
            end if
            write (summary_output_unit,'(1x,1i2,1a2,1x,1p,6d15.5)') &
                orbital_principal_n(orbital_index),orbital_symmetry_labels(orbital_index),we,wa,wb,wc,wd, &
                    orbital_occupation_weights(orbital_index)
        end do
        !
        if (selected_level_count /= 0) then
            mode = 0
            call write_scf_energies (shifted_eigenvalues,requested_level_indices,selected_level_count,mode)
            !GG         CALL ENGOUT (EVAL, IATJPO, IASPAR, ICCMIN, NCMIN, MODE)
            call write_csf_weights (.false.)
        end if
        !
        close(summary_output_unit)
        !
        return
        !
        return
        !
    end subroutine finish_summary
end module finish_summary_interface

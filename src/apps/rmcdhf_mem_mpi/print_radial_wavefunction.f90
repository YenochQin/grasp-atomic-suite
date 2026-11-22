module print_radial_wavefunction_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: print_radial_wavefunction
contains

    ! Makes a (debug) printout of wave functions. There are two modes:
    !
    !    (1) J>0  - Used as a debug option in SOLVE, wavefunctions for
    !               orbital J are printed
    !    (2) J=0  - A  printout  of the grid and all wave functions is
    !               made
    !
    ! Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992
    subroutine print_radial_wavefunction(orbital_index)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: max_radial_points
        use atomic_calculation_state, only: max_radial_points, nuclear_charge, speed_of_light_au
        use radial_grid_state, only: radial_grid
        use radial_dirac_workspace, only: max_radial_points, trial_orbital_last_grid_index, trial_large_component, &
            trial_small_component
        use orbital_csf_basis, only: orbital_count, orbital_principal_n, orbital_symmetry_labels
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        use radial_ascii_plot, only: plot_radial_orbital_ascii
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: orbital_index

        ! Local declarations
        integer :: nb2, row_count, ii, ii1, ii2, k, mfk
        real(real64) :: cbz

        cbz = speed_of_light_au/nuclear_charge
        if (orbital_index > 0) then
            ! Mode (1)
            call ensure_debug_output()
            write (angular_debug_unit, '("1",1i2,1a2," orbital:"/,2(" --------- r --------- ------- P (r) -------",&
                &" ------- Q (r) -------"))') orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
            nb2 = trial_orbital_last_grid_index/2
            if (2*nb2 == trial_orbital_last_grid_index) then
                row_count = nb2
            else
                row_count = nb2 + 1
            end if
            do ii = 1, row_count
                ii1 = ii
                ii2 = ii1 + row_count
                if (ii2 <= trial_orbital_last_grid_index) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), trial_large_component(ii1), &
                        trial_small_component(ii1), radial_grid(ii2), trial_large_component( &
                        & ii2), trial_small_component(ii2)
                else if (ii1 <= trial_orbital_last_grid_index) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), trial_large_component(ii1), &
                        trial_small_component(ii1)
                end if
            end do
            call plot_radial_orbital_ascii (trial_large_component, 1.0_real64, trial_small_component, cbz, &
                trial_orbital_last_grid_index)
        else
            ! Mode (2)
            do k = 1, orbital_count
                call ensure_debug_output()
                write (angular_debug_unit, '("1",1i2,1a2," orbital:"/,2(" --------- r --------- ------- P (r) -------",&
                    &" ------- Q (r) -------"))') orbital_principal_n(k), orbital_symmetry_labels(k)
                mfk = orbital_last_grid_index(k)
                nb2 = mfk/2
                if (2*nb2 == mfk) then
                    row_count = nb2
                else
                    row_count = nb2 + 1
                end if
                do ii = 1, row_count
                    ii1 = ii
                    ii2 = ii1 + row_count
                    if (ii2 <= mfk) then
                        call ensure_debug_output()
                        write (angular_debug_unit, '(1p,6(1x,1d21.14))') &
                            radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k), &
                            & radial_grid(ii2), radial_large_components(ii2,&
                            k), radial_small_components(ii2,k)
                    else if (ii1 <= mfk) then
                        call ensure_debug_output()
                        write (angular_debug_unit, '(1p,6(1x,1d21.14))') &
                            radial_grid(ii1), radial_large_components(ii1,k), radial_small_components(ii1,k)
                    end if
                end do
                call plot_radial_orbital_ascii (radial_large_components(:max_radial_points,k), 1.0_real64, &
                    & radial_small_components(:max_radial_points,k), cbz, orbital_last_grid_index(k))
            end do
        end if
        return
        return
    end subroutine print_radial_wavefunction
end module print_radial_wavefunction_interface

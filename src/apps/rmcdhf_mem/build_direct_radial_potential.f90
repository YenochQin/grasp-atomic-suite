! Form the direct contribution to the radial equation of the current orbital.
module build_direct_radial_potential_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_direct_radial_potential
contains

    ! This subroutine  tabulates the  potential function Y(r) (Eq (14)
    ! in  I P Grant, B J McKenzie, P H Norrington, D F Mayers, and N C
    ! Pyper, Computer  Phys  Commun  21 (1980) 211) for orbital J. The
    ! function is tabulated in the COMMON array YP.
    !
    ! Written by Farid A Parpia, at Oxford. Last update: 17 Dec 1992
    ! MPI version by Xinghong He. Last revision: 05 Aug 1998
    subroutine build_direct_radial_potential(orbital_index)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix
        use debug_output_settings, only: rmcdhf_debug_flags
        use radial_grid_state, only: radial_point_count, radial_grid
        use nuclear_potential_state, only: negative_radius_times_potential
        use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
        use scf_radial_potentials, only: orbital_y_potential
        use scf_iteration_state, only: y_potential_coefficient_count, y_potential_coefficients, y_potential_coefficient_labels
        use radial_quadrature_workspace, only: radial_workspace_b
        use hartree_yz_functions, only: evaluate_hartree_yz_functions
        use radial_ascii_plot, only: plot_radial_orbital_ascii
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index

        ! Local constants

        ! Local declarations
        integer :: index, label, k, ioy1, ioy2, nb3, row_count, ii, ii1, ii2, ii3
        real(real64) :: coeff
        ! Debug printout: composition

        if (rmcdhf_debug_flags(29) .or. rmcdhf_debug_flags(30)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/," Direct potential for ",1i2,1a2," orbital :"/,/)') &
                orbital_principal_n(orbital_index), orbital_symmetry_labels(orbital_index)
        end if
        ! Initialize array YP with the nuclear potential piece
        !
        ! Since YA() below  contains contributions from THIS node only,
        ! the initialization should be in consistence with that.

        orbital_y_potential(:radial_point_count) = negative_radius_times_potential(:radial_point_count)

        do index = 1, y_potential_coefficient_count

            ! Decode information in label
            label = y_potential_coefficient_labels(index)
            k = mod(label,orbital_index_radix)
            label = label/orbital_index_radix
            ioy1 = mod(label,orbital_index_radix)
            ioy2 = label/orbital_index_radix
            coeff = y_potential_coefficients(index)

            if (rmcdhf_debug_flags(29)) then
                call ensure_debug_output()
                write (angular_debug_unit, '(/,25x,"(",1i2,")"/,1x,1p,d21.14,"* Y    (",1i2,1a2,",",1i2,1a2,")")') &
                    k, coeff, orbital_principal_n(ioy1), orbital_symmetry_labels(ioy1), &
                    & orbital_principal_n(ioy2)&
                    , orbital_symmetry_labels(ioy2)
            end if

            call evaluate_hartree_yz_functions (k, ioy1, ioy2)                ! Accumulate contributions
            orbital_y_potential(:radial_point_count) = orbital_y_potential(:radial_point_count) - &
                coeff*radial_workspace_b(:radial_point_count)
        end do

        if (rmcdhf_debug_flags(30)) then
            call ensure_debug_output()
            write (angular_debug_unit, '(/,/,3(" --------- r --------- ------- Y (r) -------"))')
            nb3 = radial_point_count/3
            if (3*nb3 == radial_point_count) then
                row_count = nb3
            else
                row_count = nb3 + 1
            end if
            do ii = 1, row_count
                ii1 = ii
                ii2 = ii1 + row_count
                ii3 = ii2 + row_count
                if (ii3 <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') &
                        radial_grid(ii1), orbital_y_potential(ii1), radial_grid(ii2), orbital_y_potential(ii2), &
                            radial_grid(ii3), orbital_y_potential(ii3&
                        )
                else if (ii2 <= radial_point_count) then
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), orbital_y_potential(ii1), &
                        radial_grid(ii2), orbital_y_potential(ii2)
                else
                    call ensure_debug_output()
                    write (angular_debug_unit, '(1p,6(1x,1d21.14))') radial_grid(ii1), orbital_y_potential(ii1)
                end if
            end do
            call plot_radial_orbital_ascii (orbital_y_potential, 1.0_real64, orbital_y_potential, 0.0_real64, radial_point_count)
        end if
        return
        return
    end subroutine build_direct_radial_potential
end module build_direct_radial_potential_interface

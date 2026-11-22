! Add the Lagrange-multiplier terms that constrain orbital overlaps.
module apply_lagrange_constraints_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: apply_lagrange_constraints
contains

    !   This  routine  includes  the Lagrange multiplier contribution in
    !   the 'exchange' term.
    !
    !   Parameter nprocs is added so that it works for mpi and serial
    !   programs.
    !
    ! Last update: 08 Dec 1992
    ! Modified by Xinghong He. Last update: 17 Aug 1998
    subroutine apply_lagrange_constraints(orbital_index, mpi_process_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,    only: orbital_index_radix
        use atomic_calculation_state, only: working_tolerance, speed_of_light_au
        use radial_grid_state, only: radial_grid
        use orbital_lagrange_multipliers, only: lagrange_multiplier_values, lagrange_orbital_pair_labels, lagrange_multiplier_count
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source
        use scf_iteration_state, only: orbital_occupation_weights
        use radial_wavefunctions, only: orbital_last_grid_index, radial_large_components, radial_small_components
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_index
        integer, intent(in) :: mpi_process_count

        ! Local constants

        ! Local declarations
        integer :: k, iecck, l1, l2, m, mfm, i
        real(real64) :: eps, wb, wa, wari

        if (lagrange_multiplier_count == 0) return
        ! EPS = ACCY*0.1D 00
        eps = working_tolerance*0.01_real64
        ! Add contributions from off-diagonal parameters to exchange
        wb = 1.0_real64/(orbital_occupation_weights(orbital_index)*speed_of_light_au)/mpi_process_count
        do k = 1, lagrange_multiplier_count
            ! Decode index
            iecck = lagrange_orbital_pair_labels(k)
            l1 = iecck/orbital_index_radix
            l2 = iecck - orbital_index_radix*l1
            if (orbital_index == l1) then
                m = l2
            else if (orbital_index == l2) then
                m = l1
            else
                cycle
            end if
            wa = lagrange_multiplier_values(k)*wb
            if (abs(wa) < eps) cycle
            ! ADD CONTRIBUTIONS TO EXCHANGE TERMS
            mfm = orbital_last_grid_index(m)
            do i = 1, mfm
                wari = wa*radial_grid(i)
                large_exchange_source(i) = large_exchange_source(i) + wari*radial_small_components(i,m)
                small_exchange_source(i) = small_exchange_source(i) - wari*radial_large_components(i,m)
            end do
        end do
        return
    end subroutine apply_lagrange_constraints
end module apply_lagrange_constraints_interface

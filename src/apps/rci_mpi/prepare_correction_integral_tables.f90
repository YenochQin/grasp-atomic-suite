!************************************************************************
subroutine prepare_correction_integral_tables(max_orbital_two_j, atwinv)
    ! Module dependencies
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals
    use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
    use atomic_calculation_state, only: emn, max_radial_points
    use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
        radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
    use nuclear_charge_distribution, only: max_radial_points, zdist
    use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
    use coulomb_integral_table, only: frstco, vacuum_polarization_integral_indices => indoei, ncoei, ndcoea, &
        & vacuum_polarization_integral_values => valoei
    use vacuum_polarization_integral_table, only: vacuum_polarization_integral_indices => indvpi, &
        & vacuum_polarization_integral_values => valvpi, first_vacuum_polarization_lookup => frstvp, nvpi
    use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
        valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
    use kinetic_energy_integral_table, only: frstki, indkei, ndkea, nkei, valkei
    use vinti_integral_table, only: first_vinti_lookup => frstvi, indtei, ndvin, nvinti, valtei
    use orbital_csf_basis, only: csf_count, orbital_count, csf_subshell_occupations
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit

    ! Procedure interfaces
    use build_nuclear_charge_distribution_interface, only: build_nuclear_charge_distribution
    use build_vacuum_polarization_interface, only: build_vacuum_polarization
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: max_orbital_two_j
    real(real64), intent(out) :: atwinv

    ! Local declarations
    integer :: i, k
    logical :: ldbpg
    !-----------------------------------------------
    !
    frstco = .true.
    ncoei = 0

    if (ltrans) then
        !        ...Check the maximum numbers of orbtitals allowed in brint.f
        i = max_orbitals
        select case (max_orbital_two_j)
        case (11)
            i = 114
        case (12)
            i = 112
        case (13)
            i = 110
        case (14)
            i = 108
        case (15)
            i = 106
        case (16)
            i = 105
        case (17)
            i = 103
        case (18)
            i = 101
        case (19)
            i = 100
        case (21:)
            i = 90
        end select

        if (i < orbital_count) then
            write (stderr_unit, *) "In setham. The number of orbitals is too"
            write (stderr_unit, *) "large for the brint routine"
            stop
        end if

        first = .true.
        ntpi = 0
    end if
    !
    !   Initialisations for the vacuum polarisation corrections
    !
    if (lvp) then
        call build_nuclear_charge_distribution
        call build_vacuum_polarization
        zdist(2:radial_point_count) = tb(2:radial_point_count)*radial_grid_derivative(2:radial_point_count)
        first_vacuum_polarization_lookup = .true.
        nvpi = 0
    end if
    !
    !   Initialisations for nuclear translational energy corrections
    !
    if (emn > 0.0_real64) then
        atwinv = 1.0_real64/emn
        if (lnms) then
            frstki = .true.
            nkei = 0
        end if
        if (lsms) then
            first_vinti_lookup = .true.
            nvinti = 0
        end if
    else
        ! atwinv will not be used
        lnms = .false.
        lsms = .false.
    end if

    return
end subroutine prepare_correction_integral_tables

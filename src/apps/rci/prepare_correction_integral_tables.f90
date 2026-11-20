module prepare_correction_integral_tables_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: prepare_correction_integral_tables
contains

    !************************************************************************
    subroutine prepare_correction_integral_tables(max_orbital_two_j, atwinv)
        ! Module dependencies
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: max_orbitals
        use hamiltonian_correction_flags, only: include_normal_mass_shift, include_specific_mass_shift, &
            include_transverse_interaction, include_vacuum_polarization
        use atomic_calculation_state, only: nuclear_mass_in_electron_masses
        use radial_grid_state, only: radial_point_count, radial_grid_derivative
        use nuclear_charge_distribution, only: zdist
        use radial_quadrature_workspace, only: radial_workspace_b
        use coulomb_integral_table, only: first_coulomb_integral_lookup, coulomb_integral_count
        use vacuum_polarization_integral_table, only: first_vacuum_polarization_lookup, vacuum_polarization_integral_count
        use breit_integral_tables, only: breit_integral_counts, first_breit_integral_lookup
        use kinetic_energy_integral_table, only: first_kinetic_integral_lookup, kinetic_integral_count
        use vinti_integral_table, only: first_vinti_lookup, vinti_integral_count
        use orbital_csf_basis, only: orbital_count
        use standard_io_units, only: stderr_unit

        ! Procedure interfaces
        use build_nuclear_charge_distribution_interface, only: build_nuclear_charge_distribution
        use build_vacuum_polarization_interface, only: build_vacuum_polarization
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: max_orbital_two_j
        real(real64), intent(out) :: atwinv

        ! Local declarations
        integer :: i
        !-----------------------------------------------
        !
        first_coulomb_integral_lookup = .true.
        coulomb_integral_count = 0

        if (include_transverse_interaction) then
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

            first_breit_integral_lookup = .true.
            breit_integral_counts = 0
        end if
        !
        !   Initialisations for the vacuum polarisation corrections
        !
        if (include_vacuum_polarization) then
            call build_nuclear_charge_distribution
            call build_vacuum_polarization
            zdist(2:radial_point_count) = radial_workspace_b(2:radial_point_count)*radial_grid_derivative(2:radial_point_count)
            first_vacuum_polarization_lookup = .true.
            vacuum_polarization_integral_count = 0
        end if
        !
        !   Initialisations for nuclear translational energy corrections
        !
        if (nuclear_mass_in_electron_masses > 0.0_real64) then
            atwinv = 1.0_real64/nuclear_mass_in_electron_masses
            if (include_normal_mass_shift) then
                first_kinetic_integral_lookup = .true.
                kinetic_integral_count = 0
            end if
            if (include_specific_mass_shift) then
                first_vinti_lookup = .true.
                vinti_integral_count = 0
            end if
        else
            ! atwinv will not be used
            include_normal_mass_shift = .false.
            include_specific_mass_shift = .false.
        end if

        return
    end subroutine prepare_correction_integral_tables
end module prepare_correction_integral_tables_interface

module mpi_isotope_input
    implicit none (type, external)

    ! Module interface
    private
    public :: read_mpi_isotope_data
contains
    ! Read isotope data on rank 0 and broadcast nuclear charge, mass, model, spin and moments.
    subroutine read_mpi_isotope_data (isotope_file_path)

        ! Module dependencies
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_integer, mpi_rank
        use atomic_calculation_state, only: nuclear_mass_in_electron_masses, nuclear_charge
        use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
        use isotope_properties, only: nuclear_spin, nuclear_magnetic_dipole_moment, nuclear_electric_quadrupole_moment
        use isotope_data_input, only: read_isotope_data
        implicit none (type, external)

        ! Arguments
        character (len = *) :: isotope_file_path
        if (mpi_rank == 0) call read_isotope_data (isotope_file_path)

        call mpi_bcast (nuclear_charge,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (nuclear_mass_in_electron_masses,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (nuclear_model_parameters,2,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (nuclear_model_parameter_count,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (nuclear_spin,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)
        call mpi_bcast (nuclear_magnetic_dipole_moment,1,mpi_double_precision,0,mpi_comm_world, &
            mpi_error_code)
        call mpi_bcast (nuclear_electric_quadrupole_moment,1,mpi_double_precision,0,mpi_comm_world, &
            mpi_error_code)

        return
    end subroutine read_mpi_isotope_data
end module mpi_isotope_input

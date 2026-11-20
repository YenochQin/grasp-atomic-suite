module build_orbital_potential_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_orbital_potential_mpi
contains

    subroutine build_orbital_potential_mpi(eol, j, npts)
        ! Module dependencies
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use scf_radial_potentials, only: large_exchange_source, small_exchange_source, orbital_y_potential
        use mpi_state, only: mpi_allreduce, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_process_count, mpi_sum

        ! Procedure interfaces
        use assemble_orbital_potential_coefficients_interface, only: assemble_orbital_potential_coefficients
        use build_direct_radial_potential_interface, only: build_direct_radial_potential
        use build_exchange_radial_potential_interface, only: build_exchange_radial_potential
        use apply_lagrange_constraints_interface, only: apply_lagrange_constraints
        use apply_off_diagonal_integrals_interface, only: apply_off_diagonal_integrals

        ! Module dependencies

        ! Procedure interfaces
        use blas_lapack_interfaces, only: dcopy
        implicit none (type, external)

        ! Arguments
        integer  :: j
        integer  :: npts
        logical  :: eol

        ! Local declarations
        real(real64) :: tmpmpi(npts)
        !-----------------------------------------------------------------------

        call assemble_orbital_potential_coefficients (eol, j)
        call build_direct_radial_potential (j)
        call build_exchange_radial_potential (j)
        call apply_lagrange_constraints (j, mpi_process_count)
        call apply_off_diagonal_integrals

        call mpi_allreduce (orbital_y_potential, tmpmpi, npts, mpi_double_precision,  &
            mpi_sum, mpi_comm_world, mpi_error_code)
        call dcopy (npts, tmpmpi, 1, orbital_y_potential, 1)

        call mpi_allreduce (large_exchange_source, tmpmpi, npts, mpi_double_precision,  &
            mpi_sum, mpi_comm_world, mpi_error_code)
        call dcopy (npts, tmpmpi, 1, large_exchange_source, 1)

        call mpi_allreduce (small_exchange_source, tmpmpi, npts, mpi_double_precision,  &
            mpi_sum, mpi_comm_world, mpi_error_code)
        call dcopy (npts, tmpmpi, 1, small_exchange_source, 1)

        return
    end subroutine build_orbital_potential_mpi
end module build_orbital_potential_mpi_interface

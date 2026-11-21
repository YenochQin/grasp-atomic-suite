module read_al_average_options_mpi_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_al_average_options_mpi
contains

    !***********************************************************************
    !
    subroutine read_al_average_options_mpi
        ! Module dependencies
        !
        !   Interactively determines the data governing AL problem.
        !
        !   Quantities to be determined:
        !     ncmin, wt(1:ncf), ucf(1:nw), nscf, nsic, orthst
        !
        !   Call(s) to: [LIB92]: ALLOC, IQ.
        !
        !   Written by Farid A. Parpia            Last revision: 19 Dec 1992
        !   Block version by Xinghong He          Last revision: 13 Jul 1998
        !
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use array_memory_management, only: alloc

        use atomic_calculation_state, only: selected_level_count, maximum_scf_iterations, orbital_updates_per_scf_iteration, &
            level_average_weights
        use fixed_orbital_flags, only: fixed_orbital_count
        use orbital_csf_basis, only: csf_count, orbital_count
        use scf_orthogonality_control, only: enforce_orbital_orthogonality
        use scf_iteration_state, only: orbital_occupation_weights
        use standard_io_units, only: stderr_unit
        use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, mpi_rank

        ! Procedure interfaces
        use read_al_level_weights_interface, only: read_al_level_weights

        ! Module dependencies
        use csf_quantum_queries, only: csf_occupation
        implicit none (type, external)

        ! Local declarations
        integer      :: orbital_index, csf_index
        real(real64) :: sum
        !-----------------------------------------------
        !
        if (mpi_rank == 0) then
            write (stderr_unit, *) "(E)AL type calculation; H(DC) will not be ", &
                "diagonalised;"
            write (stderr_unit, *) "getald ..."
            write (stderr_unit, *) "ncf=", csf_count
        end if

        call alloc (level_average_weights, csf_count, "level_average_weights", "read_al_average_options_mpi")

        if (mpi_rank == 0) call read_al_level_weights (csf_count, level_average_weights)

        call mpi_bcast (level_average_weights(1), csf_count, mpi_double_precision, 0,        &
            mpi_comm_world, mpi_error_code)

        do orbital_index = 1, orbital_count
            sum = 0.0_real64
            do csf_index = 1, csf_count
                sum = sum + level_average_weights(csf_index)*dble(csf_occupation(orbital_index,csf_index))
            end do
            orbital_occupation_weights(orbital_index) = sum
        end do

        selected_level_count = 0
        maximum_scf_iterations = 12
        orbital_updates_per_scf_iteration = 2 + (orbital_count - fixed_orbital_count)/4
        enforce_orbital_orthogonality = .false.

        return
    end subroutine read_al_average_options_mpi
end module read_al_average_options_mpi_interface

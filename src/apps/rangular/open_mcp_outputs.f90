module open_mcp_outputs_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: open_mcp_outputs
contains

    !***********************************************************************
    subroutine open_mcp_outputs(mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
        ! Module dependencies
        use mcp_generation_options, only: mcp_output_units
        !
        ! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
        ! and MPI programs.
        !
        !   Written by Xinghong He                Last revision: 30 Jun 1998
        !
        !***********************************************************************
        use atomic_calculation_state, only: electron_count
        use csf_block_cutoffs, only: active_csf_cutoffs
        use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
        use orbital_csf_basis, only: csf_count, orbital_count

        ! Procedure interfaces
        use connect_mcp_files_interface, only: connect_mcp_files
        use read_calculation_options_interface, only: read_calculation_options
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: mpi_rank
        integer  :: mpi_process_count
        integer  :: core_subshell_count
        character  :: filehead*(*)
        character  :: block_labels(*)*8

        ! Local declarations
        integer :: k

        call connect_mcp_files (mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)

        ! DIAG, ICCUTblk, LFORDR are set in GETINF

        if (mpi_rank == 0) call read_calculation_options

        do k = 30, 32 + max_coulomb_multipole_rank
            write (mcp_output_units(k)) electron_count, csf_count, orbital_count
            write (mcp_output_units(k)) diagonal_only, active_csf_cutoffs(1), use_first_order_restriction
        end do

        do k = 30, 32 + max_coulomb_multipole_rank
            close(mcp_output_units(k))
        end do

        return
    end subroutine open_mcp_outputs
end module open_mcp_outputs_interface

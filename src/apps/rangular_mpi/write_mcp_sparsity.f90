module write_mcp_sparsity_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_mcp_sparsity
contains

    !  This routine examines lists
    !                              (IC,IR,npos)
    !  to set up the array  IENDC  required by the Davidson eigensolver
    !  of Stathopoulos and Fischer.
    !
    ! Written by Farid A. Parpia. Last revision: 10 Dec 1992
    !  Modified for block interactions by C. F. Fischer        May 1997
    !  Modified for multi-processors   by Xinghong He       03 Jul 1998
    !
    ! Currently shared by mcpblk, mcpmpi
    subroutine write_mcp_sparsity(write_coefficient_callback,nonzero_element_count,print_details,symmetry_block_index,mpi_rank, &
        & mpi_process_count,indexed_output)
        ! Module dependencies
        use mpi_indexed_mcp, only: indexed_mcp, write_indexed_sparse
        use mcp_generation_options, only: mcp_temporary_units
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc, dalloc
        use orbital_csf_basis, only: csf_count
        use standard_io_units, only: stderr_unit

        ! Procedure interfaces
        use write_mcp_sparsity_summary_mpi_interface, only: write_mcp_sparsity_summary_mpi
        implicit none (type, external)

        ! Arguments and local declarations
        type(indexed_mcp), intent(inout) :: indexed_output
        procedure(write_mcp_sparsity_summary_mpi) :: write_coefficient_callback
        integer, intent(in) :: nonzero_element_count
        integer , intent(in) :: symmetry_block_index
        integer , intent(in) :: mpi_rank
        integer , intent(in) :: mpi_process_count
        logical  :: print_details
        integer :: file_block_index, previous_column_end, previous_column_index, i, column_index, cumulative_nonzero_count
        integer, pointer :: column_end_positions(:)
        integer, pointer :: row_indices(:)
        character :: mcp_channel_label*3

        if (mpi_rank == 0) write (output_unit, *) &
            "Analysing sparse matrix array definition file ...", 30

        read (mcp_temporary_units(30)) mcp_channel_label, file_block_index
        if (symmetry_block_index /= file_block_index) then
            write (stderr_unit, *) "write_mcp_sparsity: nb = ", symmetry_block_index, ".NE. mb (=", file_block_index, ")"
            stop
        end if
        !    Allocate storage for IENDC(0:NCF)
        call alloc (column_end_positions, 0,  csf_count,"column_end_positions","write_mcp_sparsity"  )
        call alloc (row_indices, nonzero_element_count, "row_indices", "write_mcp_sparsity" )
        !   Analyse data on file 30; set up IENDC and IROW
        ! In multiprocessor environment, iendc of each node will have the
        ! same length (ncf+1); but will have its own part filled. irow is
        ! local, and its length is determined by the local parameter nnonz.

        previous_column_end = 0
        previous_column_index = 0
        do i = 1, nonzero_element_count
            read (mcp_temporary_units(30)) column_index, row_indices(i), cumulative_nonzero_count
            if (column_index /= previous_column_index) then
                column_end_positions(previous_column_index) = previous_column_end
                previous_column_index = column_index
            end if
            previous_column_end = cumulative_nonzero_count
        end do
        ! xhh - changed to suits MPI environment as well
        !       IENDC(NCF) = IEND
        !      ICLAST is zero on an empty rank; IC has no value in that case.
        column_end_positions(previous_column_index) = previous_column_end
        call write_indexed_sparse(indexed_output,symmetry_block_index, &
            column_end_positions(mpi_rank+1:csf_count:mpi_process_count),row_indices)
        !    Deallocate storage
        call dalloc (column_end_positions, "column_end_positions", "write_mcp_sparsity")
        call dalloc (row_indices, "row_indices", "write_mcp_sparsity")
        call write_coefficient_callback (print_details, nonzero_element_count, csf_count)
        return
    end subroutine write_mcp_sparsity
end module write_mcp_sparsity_interface

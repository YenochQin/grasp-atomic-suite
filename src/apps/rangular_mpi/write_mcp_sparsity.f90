module write_mcp_sparsity_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_mcp_sparsity
contains

    !***********************************************************************
    !                                                                      *
    subroutine write_mcp_sparsity(write_coefficient_callback,nonzero_element_count,print_details,symmetry_block_index,mpi_rank, &
        & mpi_process_count,output_file_prefix)
        ! Module dependencies
        use mcp_generation_options, only: mcp_temporary_units
        use, intrinsic :: iso_fortran_env, only: output_unit
        !                                                                      *
        !   This routine examines lists                                        *
        !                               (IC,IR,npos)                           *
        !   to set up the array  IENDC  required by the Davidson eigensolver   *
        !   of Stathopoulos and Fischer.                                       *
        !                                                                      *
        !   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
        !                                                                      *
        !   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
        !   Modified for block interactions by C. F. Fischer        May 1997   *
        !   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
        !
        !  Currently shared by mcpblk, mcpmpi
        !
        !***********************************************************************
        use array_memory_management, only: alloc, dalloc
        use orbital_csf_basis, only: csf_count
        use standard_io_units, only: stderr_unit

        ! Procedure interfaces
        use write_mcp_sparsity_summary_mpi_interface, only: write_mcp_sparsity_summary_mpi
        implicit none (type, external)

        ! Arguments and local declarations
        integer :: mcp_append_unit
        procedure(write_mcp_sparsity_summary_mpi) :: write_coefficient_callback
        integer, intent(in) :: nonzero_element_count
        integer , intent(in) :: symmetry_block_index
        integer , intent(in) :: mpi_rank
        integer , intent(in) :: mpi_process_count
        logical  :: print_details
        character , intent(in) :: output_file_prefix*(*)
        integer :: file_block_index, previous_column_end, previous_column_index, i, column_index, cumulative_nonzero_count, &
            & io_status
        integer, pointer :: column_end_positions(:)
        integer, pointer :: row_indices(:)
        character :: mcp_channel_label*3
        !-----------------------------------------------

        if (mpi_rank == 0) write (output_unit, *) &
            "Analysing sparse matrix array definition file ...", 30

        read (mcp_temporary_units(30)) mcp_channel_label, file_block_index
        if (symmetry_block_index /= file_block_index) then
            write (stderr_unit, *) "setsda: nb = ", symmetry_block_index, ".NE. mb (=", file_block_index, ")"
            stop
        end if
        !
        !cjb           PRINT *, ' setsda.f90 FHEAD '
        !cjb           PRINT *, ' setsda.f90 FHEAD = ', FHEAD
        !   Allocate storage for IENDC(0:NCF)
        !
        call alloc (column_end_positions, 0,  csf_count,"IENDC","SETSDA"  )
        call alloc (row_indices, nonzero_element_count, "IROW", "SETSDA" )
        !
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
        !xhh - changed to suits MPI environment as well
        !      IENDC(NCF) = IEND
        !     ICLAST is zero on an empty rank; IC has no value in that case.
        column_end_positions(previous_column_index) = previous_column_end
        !
        !   Sorting complete; rewrite to mcpXXX.30 file
        !
        !cjb           PRINT *, ' before FHEAD '
        !cjb           PRINT *, ' FHEAD = ', FHEAD
        !cjb           PRINT *, ' before OPEN '
        open(newunit=mcp_append_unit, file=output_file_prefix//".30", status="OLD", form="UNFORMATTED", iostat=io_status&
            , position="APPEND")
        !cjb           PRINT *, ' after OPEN '
        if (io_status /= 0) then
            write (stderr_unit, *) " Error when opening the file mcp.30"
            stop
        end if

        write (mcp_append_unit) "MCP", symmetry_block_index, csf_count
        write (mcp_append_unit) nonzero_element_count
        write (mcp_append_unit) (column_end_positions(i),i=mpi_rank + 1,csf_count,mpi_process_count), &
            & (row_indices(i),i=1,nonzero_element_count)
        close(mcp_append_unit)
        !cjb
        !
        !   Deallocate storage
        !
        call dalloc (column_end_positions, "IENDC", "SETSDA")
        call dalloc (row_indices, "IROW", "SETSDA")
        !
        call write_coefficient_callback (print_details, nonzero_element_count, csf_count)
        return
    end subroutine write_mcp_sparsity
end module write_mcp_sparsity_interface

! Open the MCP channels and handle the symmetry-block metadata used by the calculation.
module connect_mcp_files_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: connect_mcp_files
contains

    ! Channel 30 records the sparse Hamiltonian structure and symmetry-block metadata.
    ! Channel 31 stores one-electron T coefficients. Channels 32 and above store V(k).
    ! Each channel records the generating rank, process count, basis size, and interaction settings.
    ! Written by Farid A. Parpia. Last revision: 19 Dec 1992
    ! Modified by Xinghong He. Last revision: 06 Aug 1998
    subroutine connect_mcp_files(core_subshell_count, input_symmetry_block_count, block_labels, file_prefix)
        ! Module dependencies
        use mcp_file_metadata, only: mcp_input_units
        use array_memory_management, only: alloc
        use mpi_io_errors, only: fail_input_io_mpi
        use atomic_calculation_state, only: electron_count
        use first_order_csf_cutoffs, only: first_order_csf_cutoff
        use mcp_file_metadata, only: max_mcp_file_multipole_rank
        use rmcdhf_mcp_options, only: mcp_diagonal_only, mcp_first_order_perturbation
        use mpi_state, only: mpi_process_count, mpi_rank
        use orbital_csf_basis, only: csf_count, orbital_count
        use standard_io_units, only: stderr_unit
        use hamiltonian_symmetry_blocks, only: csf_counts_by_block, symmetry_block_count
        use integer_format, only: format_integer
        use file_open, only: open_calculation_file
        implicit none (type, external)

        ! Arguments
        integer  :: core_subshell_count
        integer , intent(in) :: input_symmetry_block_count
        character , intent(in) :: file_prefix*(*)
        character  :: block_labels(*)*8

        ! Local declarations
        integer :: sparse_input_unit
        integer :: file_prefix_length, io_status, ios, i, k, text_length, file_mpi_rank, file_mpi_process_count
        logical :: found, found1
        character :: number_text*2, mcp_channel_label*3
        character(len=len(file_prefix)+3) :: file_path

        file_prefix_length = len_trim(file_prefix)

        file_path = file_prefix(1:file_prefix_length)//".30"
        sparse_input_unit = -1
        open(newunit=sparse_input_unit, file=file_path, form="UNFORMATTED", status="OLD", iostat=io_status, &
            position="asis")
        if (io_status /= 0) call fail_input_io_mpi(io_status, sparse_input_unit, "MCP connect_mcp_files open "//trim(file_path))

        ! Check each record before using its contents.

        read (sparse_input_unit, iostat=ios) core_subshell_count, symmetry_block_count, max_mcp_file_multipole_rank
        if (ios /= 0) call fail_input_io_mpi(ios, sparse_input_unit, "MCP connect_mcp_files global header")
        if (symmetry_block_count>input_symmetry_block_count .or. symmetry_block_count<1) then
            write (stderr_unit, *) "connect_mcp_files: nblock = ", symmetry_block_count
            call fail_input_io_mpi(0, sparse_input_unit, "MCP connect_mcp_files invalid block count")
        end if

        ! cjb allocate ncfblk(0:*)
        call alloc (csf_counts_by_block, 0, symmetry_block_count , "csf_counts_by_block", "connect_mcp_files")
        csf_counts_by_block(0) = 0

        read (sparse_input_unit, iostat=ios) (csf_counts_by_block(i),i=1,symmetry_block_count)
        if (ios /= 0) call fail_input_io_mpi(ios, sparse_input_unit, "MCP connect_mcp_files block sizes")
        read (sparse_input_unit, iostat=ios) (block_labels(i),i=1,symmetry_block_count)
        if (ios /= 0) call fail_input_io_mpi(ios, sparse_input_unit, "MCP connect_mcp_files block labels")

        ! Look for other mcp files

        found = .true.
        do k = 31, 32 + max_mcp_file_multipole_rank
            call format_integer (k, number_text, text_length)
            file_path = file_prefix(1:file_prefix_length)//"."//number_text(1:2)
            inquire(file=file_path, exist=found1)
            found = found .and. found1
        end do

        if (.not.found) then
            write (stderr_unit, *) "The mcp files do not exist"
            call fail_input_io_mpi(0, sparse_input_unit, "MCP connect_mcp_files missing channel file")
        end if

        ! Open the files; check file headers

        allocate(mcp_input_units(30:32 + max_mcp_file_multipole_rank))
        mcp_input_units = -1
        mcp_input_units(30) = sparse_input_unit

        do k = 30, 32 + max_mcp_file_multipole_rank

            if (k /= 30) then
                call format_integer (k, number_text, text_length)
                file_path = file_prefix(1:file_prefix_length)//"."//number_text(1:2)
                call open_calculation_file (mcp_input_units(k), file_path, "UNFORMATTED", "OLD", io_status)
                if (io_status /= 0) call fail_input_io_mpi(io_status, mcp_input_units(k), &
                    "MCP connect_mcp_files open "//trim(file_path))
            end if

            read (mcp_input_units(k), iostat=ios) mcp_channel_label, symmetry_block_count, file_mpi_rank, file_mpi_process_count
            if (ios /= 0) call fail_input_io_mpi(ios, mcp_input_units(k), "MCP connect_mcp_files block or channel header")
            if (mpi_rank/=file_mpi_rank .or. mpi_process_count/=file_mpi_process_count) then
                write (stderr_unit, *) "mcp files were generated under different", &
                    " processor configuration."
                call fail_input_io_mpi(0, mcp_input_units(k), "MCP connect_mcp_files incompatible rank layout")
            end if

            if (mcp_channel_label /= "MCP") then
                write (stderr_unit, *) "Not a sorted GRASP92 MCP File;"
                call fail_input_io_mpi(0, mcp_input_units(k), "MCP connect_mcp_files invalid file tag")
            end if

            read (mcp_input_units(k), iostat=ios) electron_count, csf_count, orbital_count
            if (ios /= 0) call fail_input_io_mpi(ios, mcp_input_units(k), "MCP connect_mcp_files basis dimensions")
            read (mcp_input_units(k), iostat=ios) mcp_diagonal_only, first_order_csf_cutoff, mcp_first_order_perturbation
            if (ios /= 0) call fail_input_io_mpi(ios, mcp_input_units(k), "MCP connect_mcp_files channel flags")

        end do

        return
    end subroutine connect_mcp_files
end module connect_mcp_files_interface

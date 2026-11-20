!***********************************************************************
!                                                                      *
subroutine connect_mcp_files(mpi_rank, mpi_process_count, core_subshell_count, block_labels, filehead)
    ! Module dependencies
    use mcp_generation_options, only: mcp_output_units, mcp_temporary_units
    !                                                                      *
    !   Open and check the  .mcp  files. File 30 stores the structure of   *
    !   H(DC) ; file 31 stores the  T  coefficients;  files 32, 33, ...,   *
    !   store V(0), V(1), ... .                                            *
    !                                                                      *
    !   Call(s) to: [LIB92]: CONVRT, GETYN, OPENFL.                        *
    !               [GENMCP]: GETINF.                                      *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 08 Dec 1992   *
    !   MPI version by Xinghong He            Last revision: 30 Jun 1998   *
    !
    !   Used by mcpvu, mcpmpivu
    !                                                                      *
    !***********************************************************************
    use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
    use orbital_csf_basis, only: orbital_count, orbital_two_j
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use integer_format, only: convrt
    use file_open, only: openfl
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: mpi_rank
    integer, intent(in) :: mpi_process_count
    integer, intent(in) :: core_subshell_count
    character, intent(in) :: filehead*(*)
    character, intent(in) :: block_labels(*)*8

    ! Local declarations
    integer :: k, lng, text_length, open_status, i
    logical :: found, found1, getyn, yes
    character :: number_text*2
    !-----------------------------------------------
    !   Determine KMAX; this is the number of  .mcp  files for the
    !   two-electron integrals

    max_coulomb_multipole_rank = 0
    do k = 1, orbital_count
        max_coulomb_multipole_rank = max(max_coulomb_multipole_rank,orbital_two_j(k))
    end do

    !   All files  mcp.xx  are UNFORMATTED;

    allocate(mcp_output_units(30:32 + max_coulomb_multipole_rank))
    allocate(mcp_temporary_units(30:32 + max_coulomb_multipole_rank))
    lng = len_trim(filehead)
    do k = 30, 32 + max_coulomb_multipole_rank
        call convrt (k, number_text, text_length)
        call openfl (mcp_output_units(k), filehead(1:lng)//"."//number_text(1:2), "UNFORMATTED", &
            "UNKNOWN", open_status)
        if (open_status == 0) cycle
        do i = 30, k - 1
            close(mcp_output_units(i))
        end do
        write (stderr_unit, *) "Error when opening the mcp files"
        stop
    end do
    !
    !   We want to know kmax before openning other mcp files (not mcp.30)
    !   in rscf
    !
    write (mcp_output_units(30)) core_subshell_count, symmetry_block_count, max_coulomb_multipole_rank
    write (mcp_output_units(30)) (csf_counts_by_block(i),i=1,symmetry_block_count)
    write (mcp_output_units(30)) (block_labels(i),i=1,symmetry_block_count)

    do k = 30, 32 + max_coulomb_multipole_rank
        write (mcp_output_units(k)) "MCP", symmetry_block_count, mpi_rank, mpi_process_count
    end do

    return
end subroutine connect_mcp_files

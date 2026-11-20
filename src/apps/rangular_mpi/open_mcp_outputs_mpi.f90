!***********************************************************************
!                                                                      *
!cjb  myid, nprocs = NOT args
subroutine open_mcp_outputs_mpi(core_subshell_count, block_labels, filehead, descriptors)
    ! Module dependencies
    use mcp_generation_options, only: mcp_output_units
    !                                                                      *
    ! Set interaction options and write the MPI MCP file headers.          *
    !                                                                      *
    !   Written by Xinghong He                Last revision: 30 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    !cjb mpi_C
    !     USE mpi_state,    ONLY: MYID, NPROCS, ierr
    use csf_descriptor_reader, only: descriptor_source
    use mpi_csf_collectives, only: check_descriptor_status
    use csf_block_cutoffs, only: iccut
    use mcp_generation_options, only: diagonal_only, max_coulomb_multipole_rank, use_first_order_restriction
    use orbital_csf_basis, only: csf_count, orbital_count
    use atomic_calculation_state, only: electron_count
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    !cjb mpi_C
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_logical, mpi_rank

    ! Procedure interfaces
    use connect_mcp_files_interface, only: connect_mcp_files
    !     USE setmcpmpi_I
    use read_calculation_options_interface, only: read_calculation_options
    implicit none (type, external)

    ! Arguments
    !      INTEGER , INTENT(IN) :: MYID
    !      INTEGER , INTENT(IN) :: NPROCS
    integer , intent(in) :: core_subshell_count
    character , intent(in) :: filehead*(*)
    character , intent(in) :: block_labels(*)*8
    type(descriptor_source), intent(in) :: descriptors

    ! Local declarations
    !cjb mpi_C
    !cjb  INTEGER :: K, ierr
    integer :: k, input_status
    !-----------------------------------------------

    input_status = 0
    if (mpi_rank == 0) call read_calculation_options(descriptors, input_status)
    call check_descriptor_status(input_status)

    call connect_mcp_files (core_subshell_count, block_labels, filehead)

    call mpi_bcast (iccut, size(iccut), mpi_integer, 0, mpi_comm_world, mpi_error_code)
    call mpi_bcast (diagonal_only,  1,mpi_logical,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast (use_first_order_restriction,1,mpi_logical,0,mpi_comm_world,mpi_error_code)

    do k = 30, 32+max_coulomb_multipole_rank
        write (mcp_output_units(k)) electron_count,csf_count,orbital_count
        write (mcp_output_units(k)) diagonal_only,iccut(1),use_first_order_restriction
    end do

    do k = 30, 32+max_coulomb_multipole_rank
        close (mcp_output_units(k))
    end do

    return
end subroutine open_mcp_outputs_mpi

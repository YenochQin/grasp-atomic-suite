program csf_record_io
    use orbital_csf_basis, only: csf_count
    use transition_basis_blocks, only: transition_block_count
    use lodcslBio_I, only: lodcslbio
#ifdef MCP_MPI
    use mpi_state, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init
#endif
    implicit none (type, external)
    integer :: csf_input_unit
    character(32) :: header, scenario
    integer :: ncore, rank, ranks, status
    logical :: corrupt
    rank = 0
    ranks = 1
#ifdef MCP_MPI
    call MPI_Init(status)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, status)
    call MPI_Comm_size(MPI_COMM_WORLD, ranks, status)
#endif
    call get_command_argument(1, scenario)
    corrupt = scenario /= "valid" .and. rank == ranks-1
    if (corrupt) then
        open(newunit=csf_input_unit, file="truncated.c", status="old")
    else
        open(newunit=csf_input_unit, file="valid.c", status="old")
    end if
    read(csf_input_unit,"(A)") header
    if (trim(header) /= "Core subshells:") error stop "CSF fixture header"
    call lodcslbio(ncore,1,csf_input_unit)
    if (corrupt) error stop "Truncated CSF was accepted"
    if (csf_count /= 1 .or. transition_block_count /= 1 .or. ncore /= 4) error stop "Complete CSF/normal EOF"
    close(csf_input_unit)
#ifdef MCP_MPI
    call MPI_Barrier(MPI_COMM_WORLD, status)
    call MPI_Finalize(status)
#endif
    print *, "PASS CSF record input"
end program csf_record_io

program csf_record_io
    use orb_c, only: csf_count
    use blk_c, only: nblock
    use lodcslBio_I, only: lodcslbio
#ifdef MCP_MPI
    use mpi_c, only: mpi_barrier, mpi_comm_rank, mpi_comm_size, mpi_comm_world, mpi_finalize, mpi_init
#endif
    implicit none (type, external)
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
        open(21, file="truncated.c", status="old")
    else
        open(21, file="valid.c", status="old")
    end if
    read(21,"(A)") header
    if (trim(header) /= "Core subshells:") error stop "CSF fixture header"
    call lodcslbio(ncore,1)
    if (corrupt) error stop "Truncated CSF was accepted"
    if (csf_count /= 1 .or. nblock /= 1 .or. ncore /= 4) error stop "Complete CSF/normal EOF"
    close(21)
#ifdef MCP_MPI
    call MPI_Barrier(MPI_COMM_WORLD, status)
    call MPI_Finalize(status)
#endif
    print *, "PASS CSF record input"
end program csf_record_io

! Exercise the production MPI sparse-file writer without launching MPI.
program sparse_empty_rank
    use orbital_csf_basis, only: ncf => csf_count
    use setsda_I, only: setsda
    use mcp_generation_options, only: mcp_temporary_units
    implicit none (type, external)
    integer :: mcp_output_unit
    integer :: rank, block, count, columns, endpoint, row
    logical :: printing
    character(3) :: label

    allocate(mcp_temporary_units(30:30))
    ncf=1
    printing=.false.
    do rank=0,1
        ! Rank 0 owns the only column; rank 1 owns no columns or elements.
        open(newunit=mcp_output_unit, file="sparse-rank.30", status="replace", form="unformatted")
        close(mcp_output_unit)
        open(newunit=mcp_temporary_units(30), status="scratch", form="unformatted")
        write(mcp_temporary_units(30)) "MCP", 1
        if (rank == 0) write(mcp_temporary_units(30)) 1, 1, 1
        rewind(mcp_temporary_units(30))
        call setsda(check_callback,1-rank,printing,1,rank,2,"sparse-rank")
        close(mcp_temporary_units(30))

        open(newunit=mcp_output_unit, file="sparse-rank.30", status="old", form="unformatted")
        read(mcp_output_unit) label, block, columns
        if (label /= "MCP" .or. block /= 1 .or. columns /= ncf) error stop "Sparse header"
        read(mcp_output_unit) count
        if (count /= 1-rank) error stop "Sparse local count"
        if (rank == 0) then
            read(mcp_output_unit) endpoint, row
            if (endpoint /= 1 .or. row /= 1) error stop "Sparse owned column"
        else
            read(mcp_output_unit) ! Required empty record for a rank with no work.
        end if
        close(mcp_output_unit, status="delete")
    end do
    print *, "PASS sparse occupied and empty rank records"
contains
    subroutine check_callback(lprint, nnonz, n)
        logical, intent(in) :: lprint
        integer, intent(in) :: nnonz, n
        if (lprint .or. nnonz /= 1-rank .or. n /= ncf) error stop "Sparse callback"
    end subroutine check_callback
end program sparse_empty_rank

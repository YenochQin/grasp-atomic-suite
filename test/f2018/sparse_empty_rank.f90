! Exercise the production MPI sparse-file writer without launching MPI.
program sparse_empty_rank
    use orb_C, only: ncf
    use setsda_I, only: setsda
    implicit none
    integer :: rank, block, count, columns, endpoint, row
    logical :: printing
    character(3) :: label

    ncf=1
    printing=.false.
    do rank=0,1
        ! Rank 0 owns the only column; rank 1 owns no columns or elements.
        open(29, file='sparse-rank.30', status='replace', form='unformatted')
        close(29)
        open(30, status='scratch', form='unformatted')
        write(30) 'MCP', 1
        if (rank == 0) write(30) 1, 1, 1
        rewind(30)
        call setsda(check_callback,1-rank,printing,1,rank,2,'sparse-rank')
        close(30)

        open(29, file='sparse-rank.30', status='old', form='unformatted')
        read(29) label, block, columns
        if (label /= 'MCP' .or. block /= 1 .or. columns /= ncf) error stop 'Sparse header'
        read(29) count
        if (count /= 1-rank) error stop 'Sparse local count'
        if (rank == 0) then
            read(29) endpoint, row
            if (endpoint /= 1 .or. row /= 1) error stop 'Sparse owned column'
        else
            read(29) ! Required empty record for a rank with no work.
        end if
        close(29, status='delete')
    end do
    print *, 'PASS sparse occupied and empty rank records'
contains
    subroutine check_callback(lprint, nnonz, n)
        logical, intent(in) :: lprint
        integer, intent(in) :: nnonz, n
        if (lprint .or. nnonz /= 1-rank .or. n /= ncf) error stop 'Sparse callback'
    end subroutine
end program

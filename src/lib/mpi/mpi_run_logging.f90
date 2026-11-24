! Summarize MPI host and directory values on rank 0.
module mpi_run_logging
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none (type, external)

    private
    public :: report_mpi_rank_values
contains
    ! All ranks supply a string of the same declared length. Rank 0 groups equal values in first-rank order.
    subroutine report_mpi_rank_values (rank, process_count, label, local_value, as_hosts)
        use mpi, only: mpi_character, mpi_comm_world, mpi_gather

        integer, intent(in) :: rank, process_count
        character(len=*), intent(in) :: label, local_value
        logical, intent(in) :: as_hosts

        character(len=len(local_value)), allocatable :: values(:)
        integer, allocatable :: rank_groups(:), group_sizes(:)
        integer :: source, group, group_count, mpi_return_code

        allocate (values(merge(process_count, 1, rank == 0)))
        call mpi_gather (local_value, len(local_value), mpi_character, values, len(local_value), &
            mpi_character, 0, mpi_comm_world, mpi_return_code)
        if (rank /= 0) return

        allocate (rank_groups(process_count), group_sizes(process_count))
        group_sizes = 0
        do source = 1, process_count
            group = findloc (values(:source), values(source), dim=1)
            rank_groups(source) = group
            group_sizes(group) = group_sizes(group) + 1
        end do
        group_count = count(group_sizes > 0)

        if (as_hosts) then
            write (output_unit, "(2a,i0,a,i0)") label, " hosts=", group_count, ", MPI ranks=", process_count
        else if (group_count == 1) then
            write (output_unit, "(3a)") label, " ", trim(values(1))
            return
        else
            write (output_unit, "(a)") label
        end if

        do group = 1, process_count
            if (group_sizes(group) == 0) cycle
            if (as_hosts) then
                write (output_unit, "(3a,i0)") "  ", trim(values(group)), ": ranks=", group_sizes(group)
            else
                write (output_unit, "(3a)", advance="no") "  ", trim(values(group)), " (ranks "
                call write_rank_ranges (pack([(source - 1, source = 1, process_count)], rank_groups == group))
            end if
        end do
    end subroutine report_mpi_rank_values

    ! Write ascending, zero-based ranks as comma-separated single ranks and contiguous ranges.
    subroutine write_rank_ranges (ranks)
        integer, intent(in) :: ranks(:)

        integer :: first, last

        first = 1
        do last = 1, size(ranks)
            if (last < size(ranks)) then
                if (ranks(last + 1) == ranks(last) + 1) cycle
            end if

            if (first > 1) write (output_unit, "(a)", advance="no") ","
            write (output_unit, "(i0)", advance="no") ranks(first)
            if (last > first) write (output_unit, "(a,i0)", advance="no") "-", ranks(last)
            first = last + 1
        end do
        write (output_unit, "(a)") ")"
    end subroutine write_rank_ranges
end module mpi_run_logging

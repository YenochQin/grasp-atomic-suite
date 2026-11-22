! Write one startup and completion summary for an MPI application.
module mpi_run_logging
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    implicit none (type, external)

    private
    public :: begin_mpi_run, report_mpi_rank_values, report_mpi_run_finish
contains
    subroutine begin_mpi_run (rank, program_name, start_ticks)
        integer, intent(in) :: rank
        character(len=*), intent(in) :: program_name
        integer(int64), intent(out) :: start_ticks

        start_ticks = 0_int64
        if (rank /= 0) return

        call system_clock (count=start_ticks)
        write (output_unit, *)
        write (output_unit, "(a)") trim(program_name)
        call write_timestamp ("Date and Time:")
    end subroutine begin_mpi_run

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

    ! Call after collective MPI finalization so late workers are included in the elapsed time.
    subroutine report_mpi_run_finish (rank, start_ticks, program_name)
        integer, intent(in) :: rank
        integer(int64), intent(in) :: start_ticks
        character(len=*), intent(in) :: program_name

        integer(int64) :: finish_ticks, ticks_per_second, clock_limit, elapsed_ticks

        if (rank /= 0) return

        call system_clock (finish_ticks, ticks_per_second, clock_limit)
        write (output_unit, *)
        call write_timestamp ("Finish Date and Time:")
        if (ticks_per_second > 0 .and. start_ticks >= 0 .and. finish_ticks >= 0) then
            if (finish_ticks >= start_ticks) then
                elapsed_ticks = finish_ticks - start_ticks
            else
                elapsed_ticks = (clock_limit - start_ticks) + finish_ticks + 1
            end if
            write (output_unit, "(a,i0,a)") "Wall time: ", elapsed_ticks / ticks_per_second, " seconds"
        else
            write (output_unit, "(a)") "Wall time: unavailable"
        end if
        write (output_unit, "(a)") trim(program_name) // ": Execution complete."
    end subroutine report_mpi_run_finish

    subroutine write_timestamp (label)
        character(len=*), intent(in) :: label

        character(len=8) :: date_text
        character(len=10) :: time_text
        character(len=5) :: zone_text

        call date_and_time (date_text, time_text, zone_text)
        write (output_unit, "(*(a))") label, " ", date_text(:4), "-", date_text(5:6), "-", date_text(7:), &
            " ", time_text(:2), ":", time_text(3:4), ":", time_text(5:), " ", trim(zone_text)
    end subroutine write_timestamp

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

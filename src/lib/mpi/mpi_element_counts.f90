! Counts for collectives that use the default-integer MPI interface.
module mpi_element_counts
    use, intrinsic :: iso_c_binding, only: c_int
    use mpi_f08, only: mpi_count_kind
    use mpi_io_errors, only: fail_input_io_mpi
    implicit none (type, external)

    private
    public :: max_mpi_count, checked_mpi_count, mpi_gather_displacements

    integer(mpi_count_kind), parameter :: max_mpi_count = min(int(huge(0),mpi_count_kind), int(huge(0_c_int),mpi_count_kind))
contains
    integer function checked_mpi_count(count, operation) result(narrow_count)
        integer(mpi_count_kind), intent(in) :: count
        character(*), intent(in) :: operation

        if (count < 0_mpi_count_kind .or. count > max_mpi_count) &
            call fail_input_io_mpi(0, -1, operation//": element count outside default MPI range")
        narrow_count = int(count)
    end function checked_mpi_count

    subroutine mpi_gather_displacements(counts, displacements, total_count)
        integer, intent(in) :: counts(:)
        integer, intent(out) :: displacements(size(counts)), total_count

        integer :: rank_index
        integer(mpi_count_kind) :: offset

        offset = 0_mpi_count_kind
        total_count = 0
        do rank_index = 1, size(counts)
            displacements(rank_index) = total_count
            if (counts(rank_index) < 0) call fail_input_io_mpi(0, -1, "MPI gather: negative element count")
            if (int(counts(rank_index),mpi_count_kind) > max_mpi_count - offset) &
                call fail_input_io_mpi(0, -1, "MPI gather displacement: element count outside default MPI range")
            offset = offset + int(counts(rank_index), mpi_count_kind)
            total_count = int(offset)
        end do
    end subroutine mpi_gather_displacements
end module mpi_element_counts

! Payload size for allocated/associated intrinsic scalars and arrays.
! This reports eight-bit bytes, matching the suite's existing GB diagnostics.
module suite_storage
    implicit none (type, external)
    private
    public :: storage_bytes
contains
    pure function storage_bytes(values) result(bytes)
        use, intrinsic :: iso_fortran_env, only: int64
        class(*), intent(in) :: values(..)
        integer(int64) :: bytes
        bytes = storage_size(values,kind=int64)/8_int64
        if (rank(values) > 0) bytes = bytes*size(values,kind=int64)
    end function storage_bytes
end module suite_storage

! Payload size for allocated/associated intrinsic scalars and arrays.
! This reports eight-bit bytes, matching the suite's existing GB diagnostics.
module suite_storage
    use vast_kind_param, only: long
    implicit none
    private
    public :: storage_bytes
contains
    pure function storage_bytes(values) result(bytes)
        class(*), intent(in) :: values(..)
        integer(long) :: bytes
        bytes = storage_size(values,kind=long)/8_long
        if (rank(values) > 0) bytes = bytes*size(values,kind=long)
    end function storage_bytes
end module suite_storage

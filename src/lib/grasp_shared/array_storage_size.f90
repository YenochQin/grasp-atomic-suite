! Counts eight-bit payload bytes of an intrinsic scalar or array.
module array_storage_size
    implicit none (type, external)

    ! Module interface
    private
    public :: storage_bytes
contains
    pure function storage_bytes(values) result(bytes)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: int64

        ! Arguments
        class(*), intent(in) :: values(..)
        integer(int64) :: bytes
        bytes = storage_size(values,kind=int64)/8_int64
        if (rank(values) > 0) bytes = bytes*size(values,kind=int64)
    end function storage_bytes
end module array_storage_size

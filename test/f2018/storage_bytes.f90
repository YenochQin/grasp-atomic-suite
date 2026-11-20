program check_storage_bytes
    use, intrinsic :: iso_fortran_env, only: int8, int64, real32, real64
    use array_storage_size, only: storage_bytes
    implicit none (type, external)
    integer(int8) :: packed(2,3,4)
    integer :: ordinary(3)
    integer(int64) :: wide(2)
    real(real64) :: real_values(2,3), empty(0)
    if (real64 /= selected_real_kind(13) .or. real64 /= kind(1.0d0)) error stop "working precision"
    if (int64 /= selected_int_kind(18)) error stop "wide integer kind"
    if (int8 /= selected_int_kind(2)) error stop "packed byte kind"
    if (real32 /= kind(1.0)) error stop "default real kind"
    if (storage_bytes(packed) /= 24_int64) error stop "packed array size"
    if (storage_bytes(ordinary) /= 12_int64) error stop "default integer width"
    if (storage_bytes(wide) /= 16_int64) error stop "wide integer size"
    if (storage_bytes(real_values) /= 48_int64) error stop "real array size"
    if (storage_bytes(real_values(:,::2)) /= 32_int64) error stop "strided payload size"
    if (storage_bytes(empty) /= 0_int64) error stop "empty array size"
    if (storage_bytes(1_int64) /= 8_int64) error stop "scalar size"
end program check_storage_bytes

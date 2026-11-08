program check_storage_bytes
    use vast_kind_param, only: byte, long, double
    use suite_storage, only: storage_bytes
    implicit none
    integer(byte) :: packed(2,3,4)
    integer :: ordinary(3)
    integer(long) :: wide(2)
    real(double) :: real_values(2,3), empty(0)
    if (storage_bytes(packed) /= 24_long) error stop 'packed array size'
    if (storage_bytes(ordinary) /= 12_long) error stop 'default integer width'
    if (storage_bytes(wide) /= 16_long) error stop 'wide integer size'
    if (storage_bytes(real_values) /= 48_long) error stop 'real array size'
    if (storage_bytes(real_values(:,::2)) /= 32_long) error stop 'strided payload size'
    if (storage_bytes(empty) /= 0_long) error stop 'empty array size'
    if (storage_bytes(1_long) /= 8_long) error stop 'scalar size'
end program check_storage_bytes

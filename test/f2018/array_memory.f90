program array_memory
    use, intrinsic :: iso_fortran_env, only: int8, real64
    use array_memory_management, only: alloc, ralloc, dalloc
    implicit none (type, external)

    integer, pointer :: indices(:)
    integer(int8), pointer :: labels(:)
    real(real64), pointer :: values(:)

    call alloc(indices, -2, 1, "indices", "array_memory")
    indices(:) = [11, 12, 13, 14]
    call ralloc(indices, -3, 3, "indices", "array_memory")
    if (lbound(indices,1) /= -3 .or. ubound(indices,1) /= 3) error stop "grown index bounds"
    if (any(indices(-2:1) /= [11, 12, 13, 14])) error stop "grown index overlap"
    call ralloc(indices, 0, 1, "indices", "array_memory")
    if (lbound(indices,1) /= 0 .or. any(indices /= [13, 14])) error stop "shrunk index overlap"
    call ralloc(indices, 0, -1, "indices", "array_memory")
    if (size(indices) /= 0) error stop "empty indices"
    call ralloc(indices, -2, -1, "indices", "array_memory")
    if (lbound(indices,1) /= -2 .or. size(indices) /= 2) error stop "grow empty indices"
    call dalloc(indices, "indices", "array_memory")
    if (associated(indices)) error stop "index release"
    call alloc(indices, 0, 0, "indices", "array_memory")
    indices(0) = 23
    if (indices(0) /= 23) error stop "index reuse"
    call dalloc(indices, "indices", "array_memory")

    call alloc(labels, 2, 4, "labels", "array_memory")
    labels(:) = [5_int8, 6_int8, 7_int8]
    call ralloc(labels, 3, 5, "labels", "array_memory")
    if (lbound(labels,1) /= 3 .or. ubound(labels,1) /= 5) error stop "shifted label bounds"
    if (any(labels(3:4) /= [6_int8, 7_int8])) error stop "shifted label overlap"
    call dalloc(labels, "labels", "array_memory")

    call alloc(values, 0, "values", "array_memory")
    if (size(values) /= 0) error stop "empty values"
    call ralloc(values, 2, "values", "array_memory")
    values(:) = [1.25_real64, -2.5_real64]
    call ralloc(values, 4, "values", "array_memory")
    if (size(values) /= 4) error stop "grown value extent"
    if (.not. all(abs(values(:2) - [1.25_real64, -2.5_real64]) <= 0.0_real64)) error stop "grown value overlap"
    call ralloc(values, 1, "values", "array_memory")
    if (size(values) /= 1) error stop "shrunk value extent"
    if (.not. (abs(values(1) - 1.25_real64) <= 0.0_real64)) error stop "shrunk value overlap"
    call dalloc(values, "values", "array_memory")
    if (associated(values)) error stop "value release"
end program array_memory

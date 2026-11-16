program natural_orbitals_test
    use, intrinsic :: iso_fortran_env, only: real64
    use natural_orbitals
    implicit none (type, external)
    real(real64) :: matrix(23,23), occupation(23), rotation(23,23), diagonal(23,23)
    real(real64) :: p(2,23), q(2,23), origin(23), before(2,23), original_origin(23)
    character(160) :: error
    integer :: i

    ! A pure empty orbital precedes a mixed block. Upstream's ordering can
    ! overwrite this pure slot and leave another slot zero; also exercise >20.
    matrix = 0.0_real64
    do i = 4,23
        matrix(i,i) = 0.001_real64*i
    end do
    matrix(2,2) = 1.5_real64
    matrix(3,3) = 0.5_real64
    matrix(2,3) = 0.5_real64
    matrix(3,2) = 0.5_real64
    call diagonalize_density(matrix, 2, occupation, rotation, error)
    if (len_trim(error) /= 0) error stop "eigensolver failed"
    if (abs(occupation(1)) > 1e-12_real64) error stop "pure slot moved"
    if (abs(occupation(2)-(1e0_real64+sqrt(0.5_real64))) > 1e-12_real64) error stop "wrong occupation"
    diagonal = matmul(transpose(rotation),matmul(matrix,rotation))
    do i = 1,23
        diagonal(i,i) = diagonal(i,i)-occupation(i)
    end do
    if (maxval(abs(diagonal)) > 1e-12_real64) error stop "density not diagonalized"
    diagonal = matmul(transpose(rotation),rotation)
    do i = 1,23
        diagonal(i,i) = diagonal(i,i)-1e0_real64
    end do
    if (maxval(abs(diagonal)) > 1e-12_real64) error stop "rotation is not orthogonal"

    p(1,:) = [(real(i,real64),i=1,23)]
    p(2,:) = -2e0_real64*p(1,:)
    before = p
    q = 0.1_real64*p
    origin = -p(1,:)
    original_origin = origin
    call transform_radial_block(rotation, p, q, origin)
    if (any(origin < 0e0_real64)) error stop "negative origin phase"
    if (maxval(abs(p-matmul(before,rotation))) > 1e-12_real64) error stop "wrong P rotation"
    if (maxval(abs(q-0.1_real64*p)) > 1e-12_real64) error stop "inconsistent Q phase"
    if (maxval(abs(origin-matmul(original_origin,rotation))) > 1e-12_real64) error stop "inconsistent PZ phase"

    ! Repeated calls must release local workspace under -fno-automatic.
    call diagonalize_density(matrix, 1, occupation, rotation, error)
    if (len_trim(error) /= 0) error stop "dominance ordering failed"
    matrix(1:2,1:2) = reshape([1e0_real64,0.5_real64,0.5_real64,1e0_real64],[2,2])
    call diagonalize_density(matrix(1:2,1:2), 1, occupation(:2), rotation(:2,:2), error)
    if (len_trim(error) == 0) error stop "ambiguous dominance accepted"
    print *, "Natural-orbital eigensystem, ordering and phase checks passed"
end program natural_orbitals_test

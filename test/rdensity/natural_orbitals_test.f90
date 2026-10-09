program natural_orbitals_test
    use vast_kind_param, only: DOUBLE
    use natural_orbitals
    implicit none
    real(DOUBLE) :: matrix(23,23), occupation(23), rotation(23,23), diagonal(23,23)
    real(DOUBLE) :: p(2,23), q(2,23), origin(23), before(2,23), original_origin(23)
    character(160) :: error
    integer :: i

    ! A pure empty orbital precedes a mixed block. Upstream's ordering can
    ! overwrite this pure slot and leave another slot zero; also exercise >20.
    matrix = 0.0D0
    do i = 4,23
        matrix(i,i) = 0.001D0*i
    enddo
    matrix(2,2) = 1.5D0
    matrix(3,3) = 0.5D0
    matrix(2,3) = 0.5D0
    matrix(3,2) = 0.5D0
    call diagonalize_density(matrix, 2, occupation, rotation, error)
    if (len_trim(error) /= 0) error stop 'eigensolver failed'
    if (abs(occupation(1)) > 1D-12) error stop 'pure slot moved'
    if (abs(occupation(2)-(1D0+sqrt(0.5D0))) > 1D-12) error stop 'wrong occupation'
    diagonal = matmul(transpose(rotation),matmul(matrix,rotation))
    do i = 1,23
        diagonal(i,i) = diagonal(i,i)-occupation(i)
    enddo
    if (maxval(abs(diagonal)) > 1D-12) error stop 'density not diagonalized'
    diagonal = matmul(transpose(rotation),rotation)
    do i = 1,23
        diagonal(i,i) = diagonal(i,i)-1D0
    enddo
    if (maxval(abs(diagonal)) > 1D-12) error stop 'rotation is not orthogonal'

    p(1,:) = [(real(i,DOUBLE),i=1,23)]
    p(2,:) = -2D0*p(1,:)
    before = p
    q = 0.1D0*p
    origin = -p(1,:)
    original_origin = origin
    call transform_radial_block(rotation, p, q, origin)
    if (any(origin < 0D0)) error stop 'negative origin phase'
    if (maxval(abs(p-matmul(before,rotation))) > 1D-12) error stop 'wrong P rotation'
    if (maxval(abs(q-0.1D0*p)) > 1D-12) error stop 'inconsistent Q phase'
    if (maxval(abs(origin-matmul(original_origin,rotation))) > 1D-12) error stop 'inconsistent PZ phase'

    ! Repeated calls must release local workspace under -fno-automatic.
    call diagonalize_density(matrix, 1, occupation, rotation, error)
    if (len_trim(error) /= 0) error stop 'dominance ordering failed'
    matrix(1:2,1:2) = reshape([1D0,0.5D0,0.5D0,1D0],[2,2])
    call diagonalize_density(matrix(1:2,1:2), 1, occupation(:2), rotation(:2,:2), error)
    if (len_trim(error) == 0) error stop 'ambiguous dominance accepted'
    print *, 'Natural-orbital eigensystem, ordering and phase checks passed'
end program

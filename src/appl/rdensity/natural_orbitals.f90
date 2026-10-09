! Density eigensystem and radial transformation, shared by serial/MPI targets.
! Algorithm derived from GRASP rdensity/natorbnew.f90 (488e875d).
module natural_orbitals
    use vast_kind_param, only: DOUBLE
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: diagonalize_density, transform_radial_block
    interface
        recursive subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import DOUBLE
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(DOUBLE), intent(inout) :: a(lda,*), work(*)
            real(DOUBLE), intent(out) :: w(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    recursive subroutine diagonalize_density(matrix, ordering, occupation, rotation, error)
        real(DOUBLE), intent(in) :: matrix(:,:)
        integer, intent(in) :: ordering
        real(DOUBLE), intent(out) :: occupation(:), rotation(:,:)
        character(*), intent(out) :: error
        real(DOUBLE), allocatable :: values(:), vectors(:,:), work(:)
        real(DOUBLE) :: query(1)
        integer, allocatable :: order(:)
        logical, allocatable :: used(:)
        integer :: n, info, i, j, slot, lwork

        error = ''
        n = size(matrix,1)
        if (n < 1 .or. size(matrix,2) /= n .or. &
                size(occupation) /= n .or. any(shape(rotation) /= [n,n])) then
            error = 'Invalid density block dimensions'
            return
        endif
        if (.not.all(ieee_is_finite(matrix))) then
            error = 'Non-finite density matrix'
            return
        endif
        allocate(values(n), vectors(n,n), order(n), used(n))
        vectors = matrix
        call dsyev('V', 'U', n, vectors, n, values, query, -1, info)
        if (info /= 0) then
            error = 'DSYEV workspace query failed'
            return
        endif
        lwork = max(1, ceiling(query(1)))
        allocate(work(lwork))
        call dsyev('V', 'U', n, vectors, n, values, work, lwork, info)
        if (info /= 0) then
            error = 'DSYEV failed to diagonalize density matrix'
            return
        endif
        order = 0
        used = .false.
        select case(ordering)
        case(1)
            do i = 1,n
                slot = maxloc(abs(vectors(:,i)), dim=1)
                if (order(slot) /= 0) then
                    error = 'Two eigenvectors have the same dominant component; use --order=2'
                    return
                endif
                order(slot) = i
            enddo
        case(2)
            ! Preserve pure orbitals first, then fill every other slot in
            ! descending occupation order. Unlike upstream, no slot can be
            ! overwritten or left unassigned by mixed/pure eigenvectors.
            do i = 1,n
                slot = maxloc(abs(vectors(:,i)), dim=1)
                if (abs(vectors(slot,i)) < 1.0D0-1.0D-12) cycle
                order(slot) = i
                used(i) = .true.
            enddo
            j = n
            do slot = 1,n
                if (order(slot) /= 0) cycle
                do while (used(j))
                    j = j-1
                enddo
                order(slot) = j
                used(j) = .true.
            enddo
        case default
            error = 'Orbital ordering must be 1 or 2'
            return
        end select
        do i = 1,n
            occupation(i) = values(order(i))
            rotation(:,i) = vectors(:,order(i))
        enddo
    end subroutine

    recursive subroutine transform_radial_block(rotation, p, q, pz)
        real(DOUBLE), intent(inout) :: rotation(:,:), p(:,:), q(:,:), pz(:)
        integer :: i

        p = matmul(p, rotation)
        q = matmul(q, rotation)
        pz = matmul(pz, rotation)
        ! Fix the phase after ALL origin coefficients have been transformed.
        do i = 1,size(pz)
            if (pz(i) >= 0.0D0) cycle
            p(:,i) = -p(:,i)
            q(:,i) = -q(:,i)
            pz(i) = -pz(i)
            rotation(:,i) = -rotation(:,i)
        enddo
    end subroutine
end module natural_orbitals

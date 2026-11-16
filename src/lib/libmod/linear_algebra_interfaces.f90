! LP64 Fortran BLAS/LAPACK interfaces: default INTEGER and suite working precision.
! Keep these declarations consistent with the linked numerical library.
module suite_linear_algebra
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    private
    public :: dcopy, ddot, daxpy, dscal, dnrm2, dgemv, dgemm, dspevx, dswap, idamax
    public :: slamch, dlamch, lsame, xerbla
    interface
        logical function lsame(ca, cb)
            implicit none (type, external)
            character(1), intent(in) :: ca, cb
        end function lsame
        subroutine xerbla(srname, info)
            implicit none (type, external)
            character(*), intent(in) :: srname
            integer, intent(in) :: info
        end subroutine xerbla
        real(real32) function slamch(cmach)
            use, intrinsic :: iso_fortran_env, only: real32
            implicit none (type, external)
            character(1), intent(in) :: cmach
        end function slamch
        function dlamch(cmach) result(value)
            import real64
            implicit none (type, external)
            character(1), intent(in) :: cmach
            real(real64) :: value
        end function dlamch
        subroutine dswap(n, x, incx, y, incy)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx, incy
            real(real64), intent(inout) :: x(*), y(*)
        end subroutine dswap
        function idamax(n, x, incx) result(index)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx
            real(real64), intent(in) :: x(*)
            integer :: index
        end function idamax
        subroutine dcopy(n, x, incx, y, incy)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine dcopy
        function ddot(n, x, incx, y, incy) result(value)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: x(*), y(*)
            real(real64) :: value
        end function ddot
        subroutine daxpy(n, alpha, x, incx, y, incy)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx, incy
            real(real64), intent(in) :: alpha, x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine daxpy
        subroutine dscal(n, alpha, x, incx)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx
            real(real64), intent(in) :: alpha
            real(real64), intent(inout) :: x(*)
        end subroutine dscal
        function dnrm2(n, x, incx) result(value)
            import real64
            implicit none (type, external)
            integer, intent(in) :: n, incx
            real(real64), intent(in) :: x(*)
            real(real64) :: value
        end function dnrm2
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import real64
            implicit none (type, external)
            character(1), intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(real64), intent(inout) :: y(*)
        end subroutine dgemv
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import real64
            implicit none (type, external)
            character(1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(real64), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(real64), intent(inout) :: c(ldc,*)
        end subroutine dgemm
        subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info)
            import real64
            implicit none (type, external)
            character(1), intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, il, iu, ldz
            real(real64), intent(in) :: vl, vu, abstol
            real(real64), intent(inout) :: ap(*)
            integer, intent(out) :: m, info
            real(real64), intent(out) :: w(*), z(ldz,*), work(*)
            integer, intent(out) :: iwork(*), ifail(*)
        end subroutine dspevx
    end interface
end module suite_linear_algebra

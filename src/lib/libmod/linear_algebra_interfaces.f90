! LP64 Fortran BLAS/LAPACK interfaces: default INTEGER and suite DOUBLE.
! Keep these declarations consistent with the linked numerical library.
module suite_linear_algebra
    use vast_kind_param, only: double
    implicit none
    private
    public :: dcopy, ddot, daxpy, dscal, dnrm2, dgemv, dgemm, dspevx, dswap, idamax
    interface
        subroutine dswap(n, x, incx, y, incy)
            import double
            integer, intent(in) :: n, incx, incy
            real(double), intent(inout) :: x(*), y(*)
        end subroutine
        function idamax(n, x, incx) result(index)
            import double
            integer, intent(in) :: n, incx
            real(double), intent(in) :: x(*)
            integer :: index
        end function
        subroutine dcopy(n, x, incx, y, incy)
            import double
            integer, intent(in) :: n, incx, incy
            real(double), intent(in) :: x(*)
            real(double), intent(inout) :: y(*)
        end subroutine
        function ddot(n, x, incx, y, incy) result(value)
            import double
            integer, intent(in) :: n, incx, incy
            real(double), intent(in) :: x(*), y(*)
            real(double) :: value
        end function
        subroutine daxpy(n, alpha, x, incx, y, incy)
            import double
            integer, intent(in) :: n, incx, incy
            real(double), intent(in) :: alpha, x(*)
            real(double), intent(inout) :: y(*)
        end subroutine
        subroutine dscal(n, alpha, x, incx)
            import double
            integer, intent(in) :: n, incx
            real(double), intent(in) :: alpha
            real(double), intent(inout) :: x(*)
        end subroutine
        function dnrm2(n, x, incx) result(value)
            import double
            integer, intent(in) :: n, incx
            real(double), intent(in) :: x(*)
            real(double) :: value
        end function
        subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import double
            character(1), intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            real(double), intent(in) :: alpha, beta, a(lda,*), x(*)
            real(double), intent(inout) :: y(*)
        end subroutine
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import double
            character(1), intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(double), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
            real(double), intent(inout) :: c(ldc,*)
        end subroutine
        subroutine dspevx(jobz, range, uplo, n, ap, vl, vu, il, iu, abstol, m, w, z, ldz, work, iwork, ifail, info)
            import double
            character(1), intent(in) :: jobz, range, uplo
            integer, intent(in) :: n, il, iu, ldz
            real(double), intent(in) :: vl, vu, abstol
            real(double), intent(inout) :: ap(*)
            integer, intent(out) :: m, info
            real(double), intent(out) :: w(*), z(ldz,*), work(*)
            integer, intent(out) :: iwork(*), ifail(*)
        end subroutine
    end interface
end module suite_linear_algebra

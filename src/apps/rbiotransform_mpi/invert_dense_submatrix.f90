!
!
!     ------------------------------------------------------------------
!     I N V M A T
!     ------------------------------------------------------------------
!
subroutine invert_dense_submatrix(a, b, matdim, ndim)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    ! FIND INVERSE OF MATRIX A
    ! INPUT :
    !        A : MATRIX TO BE INVERTED
    !        B : SCRATCH ARRAY
    !        MATDIM : PHYSICAL DIMENSION OF MATRICES
    !        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
    !
    ! OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
    ! WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
    !
    !************************************************************************
    use, intrinsic :: iso_fortran_env, only: real64

    ! Procedure interfaces
    use invert_banded_matrix_interface, only: invert_banded_matrix
    use write_matrix_interface, only: write_matrix
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: matdim
    integer, intent(in) :: ndim
    real(real64)  :: a(matdim,matdim)
    real(real64)  :: b(matdim,matdim)

    ! Local declarations
    integer      :: itest, ntest
    real(real64) :: determ, epsil
    !-----------------------------------------------
    !
    itest = 0
    if (ndim == 1) then
        if (a(1,1) /= 0.0_real64) then
            a(1,1) = 1.0_real64/a(1,1)
        else
            itest = 1
        end if
    else
        determ = 0.0_real64
        epsil = 0.0_real64
        call invert_banded_matrix (a, b, ndim, determ, epsil, itest, matdim)
    end if
    !
    if (itest /= 0) then
        write (output_unit, "(A,I3)") " INVERSION PROBLEM NUMBER..", itest
        stop
    end if
    ntest = 0
    if (ntest /= 0) then
        write (output_unit, *) " INVERTED MATRIX "
        call write_matrix (a, ndim, ndim, matdim, matdim)
    end if
    !
    return
end subroutine invert_dense_submatrix

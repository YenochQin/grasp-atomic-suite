module invert_dense_submatrix_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: invert_dense_submatrix
contains

    !
    !
    !     ------------------------------------------------------------------
    !     I N V M A T
    !     ------------------------------------------------------------------
    !
    subroutine invert_dense_submatrix(matrix, inverse_workspace, leading_dimension, submatrix_order)
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
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use invert_banded_matrix_interface, only: invert_banded_matrix
        use write_matrix_interface, only: write_matrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: leading_dimension
        integer, intent(in) :: submatrix_order
        real(real64), intent(inout) :: matrix(leading_dimension,leading_dimension)
        real(real64)  :: inverse_workspace(leading_dimension,leading_dimension)

        ! Local declarations
        integer      :: inversion_status, trace_level
        real(real64) :: determinant, singularity_tolerance
        !-----------------------------------------------
        !
        inversion_status = 0
        if (submatrix_order == 1) then
            if (matrix(1,1) /= 0.0_real64) then
                matrix(1,1) = 1.0_real64/matrix(1,1)
            else
                inversion_status = 1
            end if
        else
            determinant = 0.0_real64
            singularity_tolerance = 0.0_real64
            call invert_banded_matrix (matrix, inverse_workspace, submatrix_order, determinant, singularity_tolerance, &
                & inversion_status, leading_dimension)
        end if
        !
        if (inversion_status /= 0) then
            write (output_unit, "(A,I3)") " INVERSION PROBLEM NUMBER..", inversion_status
            stop
        end if
        trace_level = 0
        if (trace_level /= 0) then
            write (output_unit, *) " INVERTED MATRIX "
            call write_matrix (matrix, submatrix_order, submatrix_order, leading_dimension, leading_dimension)
        end if
        !
        return
    end subroutine invert_dense_submatrix
end module invert_dense_submatrix_interface

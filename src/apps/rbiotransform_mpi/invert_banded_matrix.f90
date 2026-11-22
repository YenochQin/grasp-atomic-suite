! Invert the banded matrix used in the orbital transformation.
module invert_banded_matrix_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: invert_banded_matrix
contains

    !        DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
    !        FROM "DLYTAP".
    !
    ! *      DOUBLE PRECISION E,F
    ! *      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
    subroutine invert_banded_matrix(matrix, inverse_workspace, matrix_order, determinant, singularity_tolerance, &
                & inversion_status, leading_dimension)
        !  Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: matrix_order
        integer , intent(out) :: inversion_status
        integer , intent(in) :: leading_dimension
        real(real64) , intent(out) :: determinant
        real(real64) , intent(in) :: singularity_tolerance
        real(real64) , intent(inout) :: matrix(leading_dimension,leading_dimension)
        real(real64) , intent(inout) :: inverse_workspace(leading_dimension,leading_dimension)

        ! Local declarations
        integer :: diagonal_step, back_substitution_mode, row_index, column_index, penultimate_column, first_active_index, &
                & rotation_index, backward_column_step, backward_row_step, first_trailing_index
        real(real64) :: rotation_norm_or_value, rotation_cosine, rotation_sine, diagonal_product, largest_diagonal_pivot, &
                & smallest_diagonal_pivot, pivot_ratio_tolerance, pivot_ratio
        if (matrix_order < 2) then
            inversion_status = -1
            return
        end if
        diagonal_step = 0
        back_substitution_mode = 2
        ! SET EL = IDENTITY MATRIX
        do row_index = 1, matrix_order
            inverse_workspace(row_index,:matrix_order) = 0.0_real64
            inverse_workspace(row_index,row_index) = 1.0_real64
        end do
        ! TRIANGULARIZE A, FORM EL
        penultimate_column = matrix_order - 1
        first_active_index = 2
        do column_index = 1, penultimate_column
            do row_index = first_active_index, matrix_order
                if (matrix(row_index,column_index) == 0.0_real64) cycle
                rotation_norm_or_value = dsqrt(matrix(column_index,column_index)*matrix(column_index,column_index)+ &
                & matrix(row_index,column_index)*matrix(row_index,column_index))
                rotation_cosine = matrix(column_index,column_index)/rotation_norm_or_value
                rotation_sine = matrix(row_index,column_index)/rotation_norm_or_value
                do rotation_index = column_index, matrix_order
                    rotation_norm_or_value = rotation_cosine*matrix(column_index,rotation_index) + rotation_sine* &
                & matrix(row_index,rotation_index)
                    matrix(row_index,rotation_index) = rotation_cosine*matrix(row_index,rotation_index) - rotation_sine* &
                & matrix(column_index,rotation_index)
                    matrix(column_index,rotation_index) = rotation_norm_or_value
                end do
                do rotation_index = 1, matrix_order
                    rotation_norm_or_value = rotation_cosine*inverse_workspace(column_index,rotation_index) + rotation_sine* &
                & inverse_workspace(row_index,rotation_index)
                    inverse_workspace(row_index,rotation_index) = rotation_cosine*inverse_workspace(row_index,rotation_index) - &
                & rotation_sine*inverse_workspace(column_index,rotation_index)
                    inverse_workspace(column_index,rotation_index) = rotation_norm_or_value
                end do
            end do
            first_active_index = first_active_index + 1
        end do
        ! CALCULATE THE DETERMINANT
        diagonal_product = matrix(1,1)
        do row_index = 2, matrix_order
            diagonal_product = diagonal_product*matrix(row_index,row_index)
        end do
        determinant = diagonal_product
        ! IS MATRIX SINGULAR
        largest_diagonal_pivot = matrix(1,1)
        smallest_diagonal_pivot = matrix(1,1)
        do row_index = 2, matrix_order
            if (dabs(largest_diagonal_pivot) < dabs(matrix(row_index,row_index))) largest_diagonal_pivot = matrix(row_index, &
                & row_index)
            if (dabs(smallest_diagonal_pivot) <= dabs(matrix(row_index,row_index))) cycle
            smallest_diagonal_pivot = matrix(row_index,row_index)
        end do
        pivot_ratio_tolerance = singularity_tolerance
        if (pivot_ratio_tolerance <= 0) pivot_ratio_tolerance = 1.0e-8_real64
        pivot_ratio = smallest_diagonal_pivot/largest_diagonal_pivot
        if (abs(pivot_ratio) < pivot_ratio_tolerance) then
            inversion_status = 1
            return
        end if
        ! INVERT TRIANGULAR MATRIX
        column_index = matrix_order
        do backward_column_step = 1, matrix_order
            row_index = column_index
            diagonal_step = 1
            do backward_row_step = 1, column_index
                if (diagonal_step == 0) back_substitution_mode = 2
                if (diagonal_step == 1) then
                    back_substitution_mode = 1
                    diagonal_step = 0
                end if
                select case (back_substitution_mode)
                case default
                    matrix(row_index,column_index) = 1.0_real64/matrix(row_index,row_index)
                case (2)
                    first_trailing_index = row_index + 1
                    rotation_norm_or_value = 0.0_real64
                    rotation_norm_or_value = sum(matrix(row_index,first_trailing_index:column_index)* &
                & matrix(first_trailing_index:column_index,column_index))
                    matrix(row_index,column_index) = -rotation_norm_or_value/matrix(row_index,row_index)
                end select
                row_index = row_index - 1
            end do
            column_index = column_index - 1
        end do

        !        PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
        first_active_index = 1
        do row_index = 1, matrix_order
            do column_index = 1, matrix_order
                rotation_norm_or_value = 0.0_real64
                rotation_norm_or_value = sum(matrix(row_index,first_active_index:matrix_order)* &
                & inverse_workspace(first_active_index:matrix_order,column_index))
                inverse_workspace(row_index,column_index) = rotation_norm_or_value
            end do
            first_active_index = first_active_index + 1
        end do
        ! RECOPY EL TO A
        matrix(:matrix_order,:matrix_order) = inverse_workspace(:matrix_order,:matrix_order)
        inversion_status = 0
        return
    end subroutine invert_banded_matrix
end module invert_banded_matrix_interface

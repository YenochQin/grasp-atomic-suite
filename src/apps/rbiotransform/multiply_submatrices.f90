module multiply_submatrices_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: multiply_submatrices
contains

    !
    !     ------------------------------------------------------------------
    !     M A T M L 4
    !     ------------------------------------------------------------------
    !
    subroutine multiply_submatrices(product_matrix, left_matrix, right_matrix, product_row_count, product_column_count, &
                & left_row_count, left_column_count, right_row_count, right_column_count, &
        transpose_mode)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !
        ! MULTIPLY A AND B TO GIVE C
        !
        !     C = A * B             FOR ITRNSP = 0
        !
        !     C = A(TRANSPOSED) * B FOR ITRNSP = 1
        !
        !     C = A * B(TRANSPOSED) FOR ITRNSP = 2
        !
        !... JEPPE OLSEN, LAST REVISION JULY 24 1987
        !
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use write_matrix_interface, only: write_matrix
        use fill_vector_interface, only: fill_vector
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: product_row_count
        integer, intent(in) :: product_column_count
        integer, intent(in) :: left_row_count
        integer, intent(in) :: left_column_count
        integer, intent(in) :: right_row_count
        integer, intent(in) :: right_column_count
        integer , intent(in) :: transpose_mode
        real(real64), intent(inout) :: product_matrix(product_row_count,product_column_count)
        real(real64), intent(in) :: left_matrix(left_row_count,left_column_count)
        real(real64), intent(in) :: right_matrix(right_row_count,right_column_count)

        ! Local declarations
        integer :: trace_level, product_column, contracted_index
        real(real64) :: right_element, transposed_right_element
        !-----------------------------------------------
        !
        trace_level = 0
        if (trace_level /= 0) then
            write (output_unit, *)
            write (output_unit, *) " A AND B MATRIX FROM MATML4 "
            write (output_unit, *)
            call write_matrix (left_matrix, left_row_count, left_column_count, left_row_count, left_column_count)
            call write_matrix (right_matrix, right_row_count, right_column_count, right_row_count, right_column_count)
            write (output_unit, *) " NCROW NCCOL NAROW NACOL NBROW NBCOL "
            write (output_unit, "(6I6)") product_row_count, product_column_count, left_row_count, left_column_count, &
                & right_row_count, right_column_count
        end if
        !
        call fill_vector (product_matrix, 0.0_real64, product_row_count*product_column_count)
        !
        if (transpose_mode == 0) then
            do product_column = 1, product_column_count
                do contracted_index = 1, right_row_count
                    right_element = right_matrix(contracted_index,product_column)
                    product_matrix(:product_row_count,product_column) = product_matrix(:product_row_count,product_column) + &
                & left_matrix(:product_row_count,contracted_index)*right_element
                end do
            end do
        end if
        !
        !
        if (transpose_mode == 1) then
            !... C = A(T) * B
            do product_column = 1, product_column_count
                do contracted_index = 1, right_row_count
                    right_element = right_matrix(contracted_index,product_column)
                    product_matrix(:product_row_count,product_column) = product_matrix(:product_row_count,product_column) + &
                & left_matrix(contracted_index,:product_row_count)*right_element
                end do
            end do
        end if
        !
        if (transpose_mode == 2) then
            !... C = A*B(T)
            do product_column = 1, product_column_count
                do contracted_index = 1, right_column_count
                    transposed_right_element = right_matrix(product_column,contracted_index)
                    product_matrix(:product_row_count,product_column) = product_matrix(:product_row_count,product_column) + &
                & left_matrix(:product_row_count,contracted_index)*transposed_right_element
                end do
            end do
        end if
        !
        !
        if (trace_level /= 0) then
            write (output_unit, *)
            write (output_unit, *) " C MATRIX FROM MATML4 "
            write (output_unit, *)
            call write_matrix (product_matrix, product_row_count, product_column_count, product_row_count, product_column_count)
        end if
        !
        return
    end subroutine multiply_submatrices
end module multiply_submatrices_interface

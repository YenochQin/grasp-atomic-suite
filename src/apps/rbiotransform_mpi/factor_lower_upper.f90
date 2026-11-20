module factor_lower_upper_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: factor_lower_upper
contains

    !
    !
    !     ------------------------------------------------------------------
    !     L U L U
    !     ------------------------------------------------------------------
    !
    subroutine factor_lower_upper(input_matrix, packed_lower_factor, packed_upper_factor, matrix_order)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !
        ! LU DECOMPOSITION OF MATRIX A
        !
        !     A = L * U
        !
        ! WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
        ! UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
        !
        ! L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
        !
        !   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
        !
        !   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
        !
        ! THIS ADRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
        ! FOR L AND  OVER ROWS FOR U .
        !
        !
        ! NO PIVOTING IS DONE HERE , SO THE SCHEME GOES :
        !
        !     LOOP OVER R=1, NDIM
        !        LOOP OVER J = R, NDIM
        !          U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
        !        END OF LOOP OVER J
        !
        !        LOOP OVER I = R+1, NDIM
        !          L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
        !        END OF LOOP OVER I
        !     END OF LOOP OVER R
        !
        ! JEPPE OLSEN , OCTOBER 1988
        !
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use vector_inner_product_interface, only: vector_inner_product
        use print_matrix_symmetry_interface, only: print_matrix_symmetry
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        real(real64), intent(in) :: input_matrix(matrix_order,matrix_order)
        real(real64), intent(out) :: packed_lower_factor(*)
        real(real64), intent(inout) :: packed_upper_factor(*)

        ! Local declarations
        integer :: pivot_index, column_index, row_index, trace_level
        real(real64) :: inverse_upper_pivot
        !-----------------------------------------------
        !
        do pivot_index = 1, matrix_order
            !
            do column_index = pivot_index, matrix_order
                packed_upper_factor(column_index*(column_index-1)/2+pivot_index) = input_matrix(pivot_index,column_index) - &
                & vector_inner_product(packed_lower_factor(pivot_index*(pivot_index-1)/2+1),packed_upper_factor(column_index* &
                & (column_index-1)/2+1),pivot_index-1)
            end do
            !
            inverse_upper_pivot = 1.0_real64/packed_upper_factor(pivot_index*(pivot_index+ &
                & 1)/2)
            packed_lower_factor(pivot_index*(pivot_index+1)/2) = 1.0_real64
            do row_index = pivot_index + 1, matrix_order
                packed_lower_factor(row_index*(row_index-1)/2+pivot_index) = (input_matrix(row_index,pivot_index)- &
                & vector_inner_product(packed_lower_factor(row_index*(row_index-1)/2+1),packed_upper_factor(pivot_index* &
                & (pivot_index-1)/2+1),pivot_index-1))&
                    *inverse_upper_pivot
            end do
            !
        end do
        !
        trace_level = 0
        if (trace_level /= 0) then
            write (output_unit, *) " L MATRIX "
            call print_matrix_symmetry (packed_lower_factor, matrix_order)
            write (output_unit, *) " U MATRIX ( TRANSPOSED ) "
            call print_matrix_symmetry (packed_upper_factor, matrix_order)
        end if
        !
        return
    end subroutine factor_lower_upper
end module factor_lower_upper_interface

module combine_upper_lower_factors_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: combine_upper_lower_factors
contains

    !
    !     ------------------------------------------------------------------
    !     U L L A
    !     ------------------------------------------------------------------
    !
    subroutine combine_upper_lower_factors(input_matrix, upper_factor, lower_factor, matrix_order, workspace)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        !
        ! Obtain U L decomposition of matrix A
        !   A = U L
        !
        ! Note that L and U are returned in full matrix form
        !
        ! Quick and dirty routine, Jeppe Olsen, November 1991
        !************************************************************************
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use factor_lower_upper_interface, only: factor_lower_upper
        use fill_vector_interface, only: fill_vector
        use write_matrix_interface, only: write_matrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matrix_order
        real(real64), intent(in) :: input_matrix(matrix_order,matrix_order)
        real(real64), intent(out) :: upper_factor(matrix_order,matrix_order)
        real(real64), intent(out) :: lower_factor(matrix_order,matrix_order)
        real(real64), intent(inout) :: workspace(*)

        ! Local declarations
        integer :: next_workspace_index, reversed_matrix_start, row_index, packed_lower_start, packed_upper_start, outer_index, &
            inner_index, trace_level
        !-----------------------------------------------
        !
        !
        ! In order to change into LU form introduce the orthogonal matrix
        !    p(i,j) = delta(i,ndim-i)
        ! and rewrite
        ! A = P P A P P  = P L U P = PLP PUP
        ! where LU is an LU decomposition of PAP, since PLP is upper
        ! tringular AND PUP is lower traingular we have obtained the goal
        !
        ! 1 : PAP in scr(klPAP)
        !
        next_workspace_index = 1
        reversed_matrix_start = next_workspace_index
        next_workspace_index = next_workspace_index + matrix_order**2
        !
        do row_index = 1, matrix_order
            workspace(reversed_matrix_start-1+row_index:matrix_order*(matrix_order-1)+reversed_matrix_start-1+ &
                & row_index:matrix_order) = input_matrix(matrix_order+1-row_index,matrix_order:1:(-1))
        end do
        ! 2 : Lu decompose PAP
        packed_lower_start = next_workspace_index
        next_workspace_index = next_workspace_index + matrix_order*(matrix_order + 1)/2
        !
        packed_upper_start = next_workspace_index
        next_workspace_index = next_workspace_index + matrix_order*(matrix_order + 1)/2
        call factor_lower_upper (workspace(reversed_matrix_start), workspace(packed_lower_start), workspace(packed_upper_start), &
                & matrix_order)
        !          LULU(A,L,U,NDIM)
        ! Storage modes
        !   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
        !   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
        !
        !. 3 : Obtain U as PLP and L as PUP
        call fill_vector (upper_factor, 0.0_real64, matrix_order**2)
        call fill_vector (lower_factor, 0.0_real64, matrix_order**2)
        !
        do outer_index = 1, matrix_order
            do inner_index = 1, outer_index
                upper_factor(inner_index,outer_index) = workspace(packed_lower_start-1+(matrix_order+1-inner_index)* &
                & (matrix_order+1-inner_index-1)/2+(matrix_order+1-&
                    outer_index))
                lower_factor(outer_index,inner_index) = workspace(packed_upper_start-1+(matrix_order+1-inner_index)* &
                & (matrix_order+1-inner_index-1)/2+(matrix_order+1-&
                    outer_index))
            end do
        end do
        !
        trace_level = 0
        if (trace_level /= 0) then
            write (output_unit, *) " U and L from combine_upper_lower_factors "
            call write_matrix (upper_factor, matrix_order, matrix_order, matrix_order, matrix_order)
            call write_matrix (lower_factor, matrix_order, matrix_order, matrix_order, matrix_order)
        end if
        !

        return
    end subroutine combine_upper_lower_factors
end module combine_upper_lower_factors_interface

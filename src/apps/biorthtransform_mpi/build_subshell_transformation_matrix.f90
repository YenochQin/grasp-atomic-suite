! Build the orbital transformation for a subshell symmetry group.
module build_subshell_transformation_matrix_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: build_subshell_transformation_matrix
contains

    subroutine build_subshell_transformation_matrix(orbital_rotation, ci_transformation, workspace, orbital_count)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        ! GENERATE PER AKE'S T MATRIX FROM A
        ! ORBITAL ROTATION MATRIX X
        !
        ! T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
        ! MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
        !
        !         TL = 1 - L
        !         TU = U ** -1
        !
        ! WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
        ! X :
        !         X = L * U
        ! WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
        ! DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
        !
        ! JEPPE OLSEN OCTOBER 1988
        use, intrinsic :: iso_fortran_env, only: real64

        ! Procedure interfaces
        use factor_lower_upper_interface, only: factor_lower_upper
        use fill_vector_interface, only: fill_vector
        use write_matrix_interface, only: write_matrix
        use invert_dense_submatrix_interface, only: invert_dense_submatrix
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: orbital_count
        real(real64)  :: orbital_rotation(orbital_count,orbital_count)
        real(real64), intent(out) :: ci_transformation(orbital_count,orbital_count)
        real(real64), intent(inout) :: workspace(*)

        ! Local declarations
        integer :: trace_level, next_workspace_index, packed_lower_start, packed_upper_start, row_index, column_index
        ! DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
        trace_level = 0
        ! . Allocate local memory
        next_workspace_index = 1
        ! KLL = KFLREE
        packed_lower_start = next_workspace_index
        next_workspace_index = packed_lower_start + orbital_count*(orbital_count + 1)/2
        packed_upper_start = next_workspace_index
        next_workspace_index = packed_upper_start + orbital_count**2
        ! .LU factorize X
        call factor_lower_upper (orbital_rotation, workspace(packed_lower_start), workspace(packed_upper_start), orbital_count)
        ! .Expand U to full matrix
        call fill_vector (ci_transformation, 0.0_real64, orbital_count**2)
        do row_index = 1, orbital_count
            do column_index = row_index, orbital_count
                ci_transformation(row_index,column_index) = workspace(packed_upper_start-1+column_index*(column_index-1)/2+ &
                & row_index)
            end do
        end do
        if (trace_level >= 10) then
            write (output_unit, *) " MATRIX TO BE INVERTED "
            call write_matrix (ci_transformation, orbital_count, orbital_count, orbital_count, orbital_count)
        end if
        ! .Invert U
        call invert_dense_submatrix (ci_transformation, workspace(packed_upper_start), orbital_count, orbital_count)
        if (trace_level >= 10) then
            write (output_unit, *) " INVERTED MATRIX "
            call write_matrix (ci_transformation, orbital_count, orbital_count, orbital_count, orbital_count)
        end if
        ! .Subtract L
        do row_index = 1, orbital_count
            ci_transformation(row_index,:row_index-1) = -workspace(packed_lower_start+row_index*(row_index-1)/2:row_index-2+ &
                & packed_lower_start+row_index*(row_index-1)/2)
        end do
        if (trace_level /= 0) then
            write (output_unit, *) " INPUT X MATRIX "
            call write_matrix (orbital_rotation, orbital_count, orbital_count, orbital_count, orbital_count)
            write (output_unit, *) " T MATRIX "
            call write_matrix (ci_transformation, orbital_count, orbital_count, orbital_count, orbital_count)
        end if
        return
    end subroutine build_subshell_transformation_matrix
end module build_subshell_transformation_matrix_interface

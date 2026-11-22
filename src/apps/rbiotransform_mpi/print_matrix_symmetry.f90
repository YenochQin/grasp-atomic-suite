module print_matrix_symmetry_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: print_matrix_symmetry
contains

    subroutine print_matrix_symmetry(a, matdim)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        ! PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
        ! THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: matdim
        real(real64), intent(in) :: a(*)

        ! Local declarations
        integer :: jstart, jstop, i, j
        jstart = 1
        jstop = 0
        do i = 1, matdim
            jstart = jstart + i - 1
            jstop = jstop + i
            write (output_unit, '("0",2x,i3,5(1x,e14.7),/,(" ",5x,5(1x,e14.7)))') i, (a(j),j=jstart,jstop)
        end do
        return
        return
    end subroutine print_matrix_symmetry
end module print_matrix_symmetry_interface

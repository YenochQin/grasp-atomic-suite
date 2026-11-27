module write_matrix_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: write_matrix
contains

    subroutine write_matrix(a, nrow, ncol, nmrow, nmcol)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: nrow
        integer, intent(in) :: ncol
        integer, intent(in) :: nmrow
        integer, intent(in) :: nmcol
        real(real64), intent(in) :: a(nmrow,nmcol)

        ! Local declarations
        integer :: i, j

        do i = 1, nrow
            write (output_unit, '("0",i5,2x,4(1x,e14.7),/,(" ",7x,4(1x,e14.7)))') i, (a(i,j),j=1,ncol)
        end do

        return
    end subroutine write_matrix
end module write_matrix_interface

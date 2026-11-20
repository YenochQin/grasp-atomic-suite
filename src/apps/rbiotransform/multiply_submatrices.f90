!
!     ------------------------------------------------------------------
!     M A T M L 4
!     ------------------------------------------------------------------
!
subroutine multiply_submatrices(c, a, b, ncrow, nccol, narow, nacol, nbrow, nbcol, &
    itrnsp)
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
    integer, intent(in) :: ncrow
    integer, intent(in) :: nccol
    integer, intent(in) :: narow
    integer, intent(in) :: nacol
    integer, intent(in) :: nbrow
    integer, intent(in) :: nbcol
    integer , intent(in) :: itrnsp
    real(real64)  :: c(ncrow,nccol)
    real(real64), intent(in) :: a(narow,nacol)
    real(real64), intent(in) :: b(nbrow,nbcol)

    ! Local declarations
    integer :: ntest, j, k, i
    real(real64) :: bkj, bjk
    !-----------------------------------------------
    !
    ntest = 0
    if (ntest /= 0) then
        write (output_unit, *)
        write (output_unit, *) " A AND B MATRIX FROM MATML4 "
        write (output_unit, *)
        call write_matrix (a, narow, nacol, narow, nacol)
        call write_matrix (b, nbrow, nbcol, nbrow, nbcol)
        write (output_unit, *) " NCROW NCCOL NAROW NACOL NBROW NBCOL "
        write (output_unit, "(6I6)") ncrow, nccol, narow, nacol, nbrow, nbcol
    end if
    !
    call fill_vector (c, 0.0_real64, ncrow*nccol)
    !
    if (itrnsp == 0) then
        do j = 1, nccol
            do k = 1, nbrow
                bkj = b(k,j)
                c(:ncrow,j) = c(:ncrow,j) + a(:ncrow,k)*bkj
            end do
        end do
    end if
    !
    !
    if (itrnsp == 1) then
        !... C = A(T) * B
        do j = 1, nccol
            do k = 1, nbrow
                bkj = b(k,j)
                c(:ncrow,j) = c(:ncrow,j) + a(k,:ncrow)*bkj
            end do
        end do
    end if
    !
    if (itrnsp == 2) then
        !... C = A*B(T)
        do j = 1, nccol
            do k = 1, nbcol
                bjk = b(j,k)
                c(:ncrow,j) = c(:ncrow,j) + a(:ncrow,k)*bjk
            end do
        end do
    end if
    !
    !
    if (ntest /= 0) then
        write (output_unit, *)
        write (output_unit, *) " C MATRIX FROM MATML4 "
        write (output_unit, *)
        call write_matrix (c, ncrow, nccol, ncrow, nccol)
    end if
    !
    return
end subroutine multiply_submatrices

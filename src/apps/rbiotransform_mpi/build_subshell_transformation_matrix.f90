!
!     ------------------------------------------------------------------
!     P A M T M T
!     ------------------------------------------------------------------
!
subroutine build_subshell_transformation_matrix(x, t, work, norb)
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !
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
    !
    !-----------------------------------------------
    use, intrinsic :: iso_fortran_env, only: real64

    ! Procedure interfaces
    use factor_lower_upper_interface, only: factor_lower_upper
    use fill_vector_interface, only: fill_vector
    use write_matrix_interface, only: write_matrix
    use invert_dense_submatrix_interface, only: invert_dense_submatrix
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: norb
    real(real64)  :: x(norb,norb)
    real(real64)  :: t(norb,norb)
    real(real64), intent(inout) :: work(*)

    ! Local declarations
    integer :: ntest, klfree, kll, klu, i, j
    !-----------------------------------------------
    ! DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
    !
    ntest = 0
    !. Allocate local memory
    klfree = 1
    !     KLL = KFLREE
    kll = klfree
    klfree = kll + norb*(norb + 1)/2
    klu = klfree
    klfree = klu + norb**2
    !.LU factorize X
    call factor_lower_upper (x, work(kll), work(klu), norb)
    !.Expand U to full matrix
    call fill_vector (t, 0.0_real64, norb**2)
    do i = 1, norb
        do j = i, norb
            t(i,j) = work(klu-1+j*(j-1)/2+i)
        end do
    end do
    if (ntest >= 10) then
        write (output_unit, *) " MATRIX TO BE INVERTED "
        call write_matrix (t, norb, norb, norb, norb)
    end if
    !.Invert U
    call invert_dense_submatrix (t, work(klu), norb, norb)
    if (ntest >= 10) then
        write (output_unit, *) " INVERTED MATRIX "
        call write_matrix (t, norb, norb, norb, norb)
    end if
    !.Subtract L
    do i = 1, norb
        t(i,:i-1) = -work(kll+i*(i-1)/2:i-2+kll+i*(i-1)/2)
    end do
    !
    if (ntest /= 0) then
        write (output_unit, *) " INPUT X MATRIX "
        call write_matrix (x, norb, norb, norb, norb)
        write (output_unit, *) " T MATRIX "
        call write_matrix (t, norb, norb, norb, norb)
    end if
    !
    return
end subroutine build_subshell_transformation_matrix

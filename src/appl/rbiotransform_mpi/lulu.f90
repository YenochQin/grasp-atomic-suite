!
!
!     ------------------------------------------------------------------
!     L U L U
!     ------------------------------------------------------------------
!
      subroutine lulu(a, l, u, ndim)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use inprod_i
      use prsym_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ndim
      real(real64), intent(in) :: a(ndim,ndim)
      real(real64)  :: l(*)
      real(real64)  :: u(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: r, j, i, ntest
      real(real64) :: xfaci
!-----------------------------------------------
!
      do r = 1, ndim
!
         do j = r, ndim
            u(j*(j-1)/2+r) = a(r,j) - inprod(l(r*(r-1)/2+1),u(j*(j-1)/2+1),r-1)
         end do
!
         xfaci = 1.0_real64/u(r*(r+1)/2)
         l(r*(r+1)/2) = 1.0_real64
         do i = r + 1, ndim
            l(i*(i-1)/2+r) = (a(i,r)-inprod(l(i*(i-1)/2+1),u(r*(r-1)/2+1),r-1))&
               *xfaci
         end do
!
      end do
!
      ntest = 0
      if (ntest /= 0) then
         write (6, *) ' L MATRIX '
         call prsym (l, ndim)
         write (6, *) ' U MATRIX ( TRANSPOSED ) '
         call prsym (u, ndim)
      endif
!
      return
      end subroutine lulu

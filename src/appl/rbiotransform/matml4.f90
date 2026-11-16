!
!     ------------------------------------------------------------------
!     M A T M L 4
!     ------------------------------------------------------------------
!
      subroutine matml4(c, a, b, ncrow, nccol, narow, nacol, nbrow, nbcol, &
         itrnsp)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use wrtmat_i, only: wrtmat
      use setvec_i, only: setvec
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntest, j, k, i
      real(real64) :: bkj, bjk
!-----------------------------------------------
!
      ntest = 0
      if (ntest /= 0) then
         write (6, *)
         write (6, *) " A AND B MATRIX FROM MATML4 "
         write (6, *)
         call wrtmat (a, narow, nacol, narow, nacol)
         call wrtmat (b, nbrow, nbcol, nbrow, nbcol)
         write (6, *) " NCROW NCCOL NAROW NACOL NBROW NBCOL "
         write (6, "(6I6)") ncrow, nccol, narow, nacol, nbrow, nbcol
      end if
!
      call setvec (c, 0.0_real64, ncrow*nccol)
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
         write (6, *)
         write (6, *) " C MATRIX FROM MATML4 "
         write (6, *)
         call wrtmat (c, ncrow, nccol, ncrow, nccol)
      end if
!
      return
      end subroutine matml4

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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use wrtmat_i
      use setvec_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncrow
      integer  :: nccol
      integer  :: narow
      integer  :: nacol
      integer  :: nbrow
      integer  :: nbcol
      integer , intent(in) :: itrnsp
      real(double)  :: c(ncrow,nccol)
      real(double)  :: a(narow,nacol)
      real(double)  :: b(nbrow,nbcol)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntest, j, k, i
      real(double) :: bkj, bjk
!-----------------------------------------------
!
      ntest = 0
      if (ntest /= 0) then
         write (6, *)
         write (6, *) ' A AND B MATRIX FROM MATML4 '
         write (6, *)
         call wrtmat (a, narow, nacol, narow, nacol)
         call wrtmat (b, nbrow, nbcol, nbrow, nbcol)
         write (6, *) ' NCROW NCCOL NAROW NACOL NBROW NBCOL '
         write (6, '(6I6)') ncrow, nccol, narow, nacol, nbrow, nbcol
      endif
!
      call setvec (c, 0.0d0, ncrow*nccol)
!
      if (itrnsp == 0) then
         do j = 1, nccol
            do k = 1, nbrow
               bkj = b(k,j)
               c(:ncrow,j) = c(:ncrow,j) + a(:ncrow,k)*bkj
            end do
         end do
      endif
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
      endif
!
      if (itrnsp == 2) then
!... C = A*B(T)
         do j = 1, nccol
            do k = 1, nbcol
               bjk = b(j,k)
               c(:ncrow,j) = c(:ncrow,j) + a(:ncrow,k)*bjk
            end do
         end do
      endif
!
!
      if (ntest /= 0) then
         write (6, *)
         write (6, *) ' C MATRIX FROM MATML4 '
         write (6, *)
         call wrtmat (c, ncrow, nccol, ncrow, nccol)
      endif
!
      return
      end subroutine matml4

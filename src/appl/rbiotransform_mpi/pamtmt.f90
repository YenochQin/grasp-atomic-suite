!
!     ------------------------------------------------------------------
!     P A M T M T
!     ------------------------------------------------------------------
!
      subroutine pamtmt(x, t, work, norb)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lulu_i, only: lulu
      use setvec_i, only: setvec
      use wrtmat_i, only: wrtmat
      use invmat_i, only: invmat
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: norb
      real(real64)  :: x(norb,norb)
      real(real64)  :: t(norb,norb)
      real(real64), intent(inout) :: work(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
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
      call lulu (x, work(kll), work(klu), norb)
!.Expand U to full matrix
      call setvec (t, 0.0_real64, norb**2)
      do i = 1, norb
         do j = i, norb
            t(i,j) = work(klu-1+j*(j-1)/2+i)
         end do
      end do
      if (ntest >= 10) then
         write (6, *) " MATRIX TO BE INVERTED "
         call wrtmat (t, norb, norb, norb, norb)
      end if
!.Invert U
      call invmat (t, work(klu), norb, norb)
      if (ntest >= 10) then
         write (6, *) " INVERTED MATRIX "
         call wrtmat (t, norb, norb, norb, norb)
      end if
!.Subtract L
      do i = 1, norb
         t(i,:i-1) = -work(kll+i*(i-1)/2:i-2+kll+i*(i-1)/2)
      end do
!
      if (ntest /= 0) then
         write (6, *) " INPUT X MATRIX "
         call wrtmat (x, norb, norb, norb, norb)
         write (6, *) " T MATRIX "
         call wrtmat (t, norb, norb, norb, norb)
      end if
!
      return
      end subroutine pamtmt

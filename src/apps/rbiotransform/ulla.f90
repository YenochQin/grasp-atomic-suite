!
!     ------------------------------------------------------------------
!     U L L A
!     ------------------------------------------------------------------
!
      subroutine ulla(a, u, l, ndim, scr)
      use, intrinsic :: iso_fortran_env, only: output_unit
!
! Obtain U L decomposition of matrix A
!   A = U L
!
! Note that L and U are returned in full matrix form
!
! Quick and dirty routine, Jeppe Olsen, November 1991
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use lulu_i, only: lulu
      use setvec_i, only: setvec
      use wrtmat_i, only: wrtmat
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ndim
      real(real64), intent(in) :: a(ndim,ndim)
      real(real64)  :: u(ndim,ndim)
      real(real64)  :: l(ndim,ndim)
      real(real64)  :: scr(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: klfree, klpap, i, j, kll, klu, imax, imin, ntest
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
      klfree = 1
      klpap = klfree
      klfree = klfree + ndim**2
!
      do i = 1, ndim
         scr(klpap-1+i:ndim*(ndim-1)+klpap-1+i:ndim) = a(ndim+1-i,ndim:1:(-1))
      end do
! 2 : Lu decompose PAP
      kll = klfree
      klfree = klfree + ndim*(ndim + 1)/2
!
      klu = klfree
      klfree = klfree + ndim*(ndim + 1)/2
      call lulu (scr(klpap), scr(kll), scr(klu), ndim)
!          LULU(A,L,U,NDIM)
! Storage modes
!   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
!   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
!
!. 3 : Obtain U as PLP and L as PUP
      call setvec (u, 0.0_real64, ndim**2)
      call setvec (l, 0.0_real64, ndim**2)
!
      do imax = 1, ndim
         do imin = 1, imax
            u(imin,imax) = scr(kll-1+(ndim+1-imin)*(ndim+1-imin-1)/2+(ndim+1-&
               imax))
            l(imax,imin) = scr(klu-1+(ndim+1-imin)*(ndim+1-imin-1)/2+(ndim+1-&
               imax))
         end do
      end do
!
      ntest = 0
      if (ntest /= 0) then
         write (output_unit, *) " U and L from Ulla "
         call wrtmat (u, ndim, ndim, ndim, ndim)
         call wrtmat (l, ndim, ndim, ndim, ndim)
      end if
!

      return
      end subroutine ulla

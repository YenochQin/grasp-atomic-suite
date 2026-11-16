!
!
!     ------------------------------------------------------------------
!     I N V M A T
!     ------------------------------------------------------------------
!
      subroutine invmat(a, b, matdim, ndim)
! FIND INVERSE OF MATRIX A
! INPUT :
!        A : MATRIX TO BE INVERTED
!        B : SCRATCH ARRAY
!        MATDIM : PHYSICAL DIMENSION OF MATRICES
!        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
!
! OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
! WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use bndinv_i
      use wrtmat_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: matdim
      integer, intent(in) :: ndim
      real(real64)  :: a(matdim,matdim)
      real(real64)  :: b(matdim,matdim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: itest, ntest
      real(real64) :: determ, epsil
!-----------------------------------------------
!
      itest = 0
      if (ndim == 1) then
         if (a(1,1) /= 0.0_real64) then
            a(1,1) = 1.0_real64/a(1,1)
         else
            itest = 1
         endif
      else
         determ = 0.0_real64
         epsil = 0.0_real64
         call bndinv (a, b, ndim, determ, epsil, itest, matdim)
      endif
!
      if (itest /= 0) then
         write (6, '(A,I3)') ' INVERSION PROBLEM NUMBER..', itest
         stop
      endif
      ntest = 0
      if (ntest /= 0) then
         write (6, *) ' INVERTED MATRIX '
         call wrtmat (a, ndim, ndim, matdim, matdim)
      endif
!
      return
      end subroutine invmat

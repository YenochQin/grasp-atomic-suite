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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use bndinv_i
      use wrtmat_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: matdim
      integer  :: ndim
      real(double)  :: a(matdim,matdim)
      real(double)  :: b(matdim,matdim)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: itest, ntest
      real(double) :: determ, epsil
!-----------------------------------------------
!
      itest = 0
      if (ndim == 1) then
         if (a(1,1) /= 0.0d0) then
            a(1,1) = 1.0d0/a(1,1)
         else
            itest = 1
         endif
      else
         determ = 0.0d0
         epsil = 0.0d0
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

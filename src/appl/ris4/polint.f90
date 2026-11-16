!************************************************************************
!
subroutine polint(xa,ya,dens)
!                                                                             *
!   This routine uses interpolating polynomial d(r) = d0 + d2*r^2 + d4*r^4    *
!   to extrapolate electron density at r = 0.                                 *
!                                                                             *
!   Written by Jorgen Ekman Jul. 2016                                         *
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none
      real(double), intent(out) :: dens
      real(double), dimension(3),  intent(in) :: xa, ya
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(3,3) :: m, mi
      real(double), dimension(3) :: rm, pm

      real(double) :: mdet
!-----------------------------------------------

      rm(1) = ya(1)
      rm(2) = ya(2)
      rm(3) = ya(3)

      m(1,1) = 1.0d0
      m(1,2) = xa(1)**2.0d0
      m(1,3) = xa(1)**4.0d0
      m(2,1) = 1.0d0
      m(2,2) = xa(2)**2.0d0
      m(2,3) = xa(2)**4.0d0
      m(3,1) = 1.d0
      m(3,2) = xa(3)**2.0d0
      m(3,3) = xa(3)**4.0d0
      mdet = m(1,1)*m(2,2)*m(3,3)-m(1,1)*m(2,3)*m(3,2)- &
           m(1,2)*m(2,1)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+ &
           m(1,3)*m(2,1)*m(3,2)-m(1,3)*m(2,2)*m(3,1)

!     DETERMINE INVERSE OF MATRIX
      mi(1,1) = (m(2,2)*m(3,3)-m(2,3)*m(3,2))/mdet
      mi(1,2) = (m(1,3)*m(3,2)-m(1,2)*m(3,3))/mdet
      mi(1,3) = (m(1,2)*m(2,3)-m(1,3)*m(2,2))/mdet

      mi(2,1) = (m(2,3)*m(3,1)-m(2,1)*m(3,3))/mdet
      mi(2,2) = (m(1,1)*m(3,3)-m(1,3)*m(3,1))/mdet
      mi(2,3) = (m(1,3)*m(2,1)-m(1,1)*m(2,3))/mdet

      mi(3,1) = (m(2,1)*m(3,2)-m(2,2)*m(3,1))/mdet
      mi(3,2) = (m(1,2)*m(3,1)-m(1,1)*m(3,2))/mdet
      mi(3,3) = (m(1,1)*m(2,2)-m(1,2)*m(2,1))/mdet

!     DETERMINE PARAMETERS FROM FIT
      pm(1) = mi(1,1)*rm(1)+mi(1,2)*rm(2)+mi(1,3)*rm(3)
      pm(2) = mi(2,1)*rm(1)+mi(2,2)*rm(2)+mi(2,3)*rm(3)
      pm(3) = mi(3,1)*rm(1)+mi(3,2)*rm(2)+mi(3,3)*rm(3)

!     Finally RHO(0) in au^{-3} is defined
      dens = pm(1)

      return
    end subroutine polint

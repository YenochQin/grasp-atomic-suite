module coulom_i
    implicit none (type, external)
    private
    public :: coulom
contains
!*******************************************************************
!                                                                  *
      subroutine coulom(j1,j2,j3,j4,l1,l2,l3,l4,k,aa)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 01  ------------   *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
!                                                                  *
!                          k   k+1  (k) (k)                        *
!     (n l j T  n l j T ::r  / r  ( C   C )::n l j T  n l j T )    *
!       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
!                                                                  *
!     SUBROUTINE CALLED:  CRE                                      *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only:  eps, zero
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      use cre_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)        :: j1,j2,j3,j4,l1,l2,l3,l4,k
      real(real64), intent(out)  :: aa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ifaz
!-----------------------------------------------
      aa=zero
      if(ittk(l1,l3,k) == 0)return
      if(ittk(l2,l4,k) == 0)return
      i=(2*k+1)/2
      aa=cre (j1,i,j3)
      if(dabs(aa) < eps)return
      aa=aa*cre (j2,i,j4)
      if(dabs(aa) < eps)return
      ifaz=l3-2*k-l1+l4-l2
      if((ifaz/4)*4 /= ifaz)aa=-aa
      return
      end subroutine coulom
end module coulom_i

!*******************************************************************
!                                                                  *
      integer function izas1(ib,qb,ik,qk)
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             March    1995  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use cons_c,          only: tenth, two
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in) :: ib, ik
      real(double), intent(in) :: qb, qk
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iqb, iqk
!-----------------------------------------------
      izas1=0
      iqb=two*dabs(qb)+tenth
      if(iqb > ib)return
      if(mod(ib+iqb,2) /= 0)return
      iqk=two*dabs(qk)+tenth
      if(iqk > ik)return
      if(mod(ik+iqk,2) /= 0)return
      izas1=1
      return
      end  function izas1

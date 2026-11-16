!*******************************************************************
!                                                                  *
      subroutine diaga1(ja1,ka,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 01  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
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
      use cons_c,          only: zero
      use m_c,             only: jlist, jjq1, jjq2, jjc1, jjc2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja1, ka, ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(real64), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij1,ia1,ib1,k1,j1,it1,it1s,ifaz,ll1
      real(real64) :: a1
!-----------------------------------------------
      recc = zero
      ij1=jlist(ja1)
      ia1=jjq1(3,ij1)-1
      ib1=jjq2(3,ij1)-1
      if(ja1 <= 2) then
        ll1=jlist(1)
        if(ja1 == 1)ll1=jlist(2)
        j1=jjq1(3,ll1)-1
        it1=jjc1(1)-1
        it1s=jjc2(1)-1
      else
        k1=ja1-2
        j1=jjc1(k1)-1
        k1=ja1-1
        it1=jjc1(k1)-1
        it1s=jjc2(k1)-1
      end if
      if(ire /= 0) then
        call sixj(ka,ib1,ia1,j1,it1,it1s,0,a1)
        a1=a1*dsqrt(dble((ia1+1)*(it1s+1)))
        ifaz=j1+it1+ib1+ka
        if((ifaz/4)*4 /= ifaz)a1=-a1
        recc=a1
        iat=1
        if(ja1 /= 1)return
        ifaz=ia1+ib1+2*j1-it1-it1s
        if((ifaz/4)*4 /= ifaz)recc=-recc
      else
        if(ixjtik(ka,ib1,ia1,j1,it1,it1s) == 0)return
        iat=1
      end if
      return
      end subroutine diaga1

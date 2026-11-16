!*******************************************************************
!                                                                  *
      subroutine diaga2(ja1,ja2,ka,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 02  --------------  *
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
      use vast_kind_param, only: double
      use cons_c,          only: zero
      use m_c,             only: jlist, jjq1, jjq2, jjc1, jjc2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use sixj_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja1, ja2, ka, ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(double), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij1,ij2,ia1,ia2,ib1,ib2,n1,n2,j2,it2,it2s,ifaz
      real(double) :: a2
!-----------------------------------------------
      recc = zero
      ij1=jlist(ja1)
      ij2=jlist(ja2)
      ia1=jjq1(3,ij1)-1
      ia2=jjq1(3,ij2)-1
      ib1=jjq2(3,ij1)-1
      ib2=jjq2(3,ij2)-1
      if(ja1 == 1.and.ja2 == 2) then
        it2=ia1
        it2s=ib1
        j2=jjc1(1)-1
      else
        n1=ja2-1
        j2=jjc1(n1)-1
        n2=ja2-2
        it2=jjc1(n2)-1
        it2s=jjc2(n2)-1
      end if
      if(ire /= 0) then
        call sixj(ka,ib2,ia2,j2,it2,it2s,0,a2)
        recc=a2*dsqrt(dble((ib2+1)*(it2+1)))
        ifaz=j2+it2s+ia2+ka
        if((ifaz/4)*4 /= ifaz)recc=-recc
        iat=1
      else
        if(ixjtik(ka,ib2,ia2,j2,it2,it2s) == 0)return
        iat=1
      end if
      return
      end subroutine diaga2

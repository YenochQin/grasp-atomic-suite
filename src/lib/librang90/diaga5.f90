!*******************************************************************
!                                                                  *
      subroutine diaga5(npeelgg,ja1,ka,ire,iat,recc)
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
      use vast_kind_param, only: double
      use cons_c,          only: zero
      use m_c,             only: jlist, jjc1, jjc2, jjq1, jjq2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use sixj_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: npeelgg,ja1,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(double), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: iti1,iti1s,ij1,iti,itis,ji
      real(double) :: a3
!-----------------------------------------------
      recc = zero
      iti1=jjc1(npeelgg-1)-1
      iti1s=jjc2(npeelgg-1)-1
      ij1=jlist(npeelgg)
      if(ja1 == npeelgg) then
        iti=jjq1(3,ij1)-1
        itis=jjq2(3,ij1)-1
        ji=jjc1(npeelgg-2)-1
      else
        ji=jjq1(3,ij1)-1
        iti=jjc1(npeelgg-2)-1
        itis=jjc2(npeelgg-2)-1
      end if
      if(ire == 0) then
        if(ixjtik(ka,itis,iti,ji,iti1,iti1s) /= 0) iat=1
      else
        call sixj(ka,itis,iti,ji,iti1,iti1s,0,a3)
        recc=a3*dsqrt(dble((iti+1)*(iti1s+1)))
        if(mod(ka+ji+itis+iti1,4) /= 0)recc=-recc
        iat=1
        if(ja1 == npeelgg)return
        if(mod(iti+itis-iti1s-iti1+2*ji,4) /= 0)recc=-recc
      end if
      return
      end subroutine diaga5

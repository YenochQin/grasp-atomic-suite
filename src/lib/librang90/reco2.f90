!*******************************************************************
!                                                                  *
      subroutine reco2(ja1,ja2,ka,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: one
      use m_c,             only: jlist, jjq1, jjq2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use diaga1_i
      use diaga2_i
      use diaga3_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja1,ja2,ka,ire
      integer, intent(out)      :: iat
      real(real64), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ia1,ia2,ib1,ib2,ij1,ij2,iskr
      real(real64) :: s, ss, re
!-----------------------------------------------
      iat=0
      ij1=jlist(ja1)
      ij2=jlist(ja2)
      s=dble(jjq1(3,ij1))
      ss=dble(jjq1(3,ij2))
      ss=s*ss
      recc=one/dsqrt(ss)
      if(ire == 0) then
        iat=0
      else if(ka /= 0) then
        iat=0
      else
        iat=1
        return
      end if
      ia1=jjq1(3,ij1)-1
      ia2=jjq1(3,ij2)-1
      ib1=jjq2(3,ij1)-1
      ib2=jjq2(3,ij2)-1
!
      call diaga2(ja1,ja2,ka,ire,iat,re)
      if(iat == 0)return
      recc=re*recc*dsqrt(dble(ia2+1))/dsqrt(dble((ka+1)*(ib2+1)))
      if(ja1 == 1.and.ja2 == 2)return
!
      iat=0
      call diaga1(ja1,ka,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      iskr=ja2-ja1
      if(ja1 == 1)iskr=ja2-1-ja1
      if(iskr <= 1)return
!
      iat=0
      call diaga3(ja1,ja2,ka,ire,iat,re)
      recc=re*recc
      return
      end subroutine reco2

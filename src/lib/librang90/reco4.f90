!*******************************************************************
!                                                                  *
      subroutine reco4(ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 09  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
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
      use diaga4_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ja1,ja2,ja3,ja4,k1,k2,k3,k4,ka,ire
      integer, intent(out) :: iat
      real(real64), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ia4, ib4, ij1, ij2, ij3, ij4, iskr
      real(real64) :: s, s1, s2, s3, s4, re
!-----------------------------------------------
      ij1=jlist(ja1)
      ij2=jlist(ja2)
      ij3=jlist(ja3)
      ij4=jlist(ja4)
      s1=jjq1(3,ij1)
      s2=jjq1(3,ij2)
      s3=jjq1(3,ij3)
      s4=jjq1(3,ij4)
      s=s1*s2*s3*s4
      recc=one/dsqrt(s)
      ia4=jjq1(3,ij4)-1
      ib4=jjq2(3,ij4)-1
      recc=recc*dsqrt(dble(ia4+1))/dsqrt(dble((k4+1)*(ib4+1)))
!
      iskr=ja3-ja2
      if(iskr > 1) then
        iat=0
        call diaga3(ja2,ja3,ka,ire,iat,re)
        if(iat == 0)return
        recc=re*recc
      end if
!
      iskr=ja4-ja3
      if(iskr > 1) then
        iat=0
        call diaga3(ja3,ja4,k4,ire,iat,re)
        if(iat == 0)return
        recc=re*recc
      end if
!
      iat=0
      call diaga2(ja1,ja4,k4,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
!
      iat=0
      call diaga4(ja1,ja2,k1,k2,ka,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
!
      iat=0
      call diaga4(ja2,ja3,ka,k3,k4,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      if(ja1 == 1.and.ja2 == 2)return
!
      iat=0
      call diaga1(ja1,k1,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
!
      iskr=ja2-ja1
      if(ja1 == 1)iskr=ja2-1-ja1
      if(iskr <= 1)return
      iat=0
      call diaga3(ja1,ja2,k1,ire,iat,re)
      recc=re*recc
      return
      end subroutine reco4

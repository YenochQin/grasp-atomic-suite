!*******************************************************************
!                                                                  *
      subroutine recop2(ns,ja1,ja2,k1,k2,ka,ire,iat,recc)
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
      use vast_kind_param, only: double
      use cons_c,          only: one
      use m_c,             only: jlist, jjq1, jjq2, npeel
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use diaga1_i
      use diaga3_i
      use diaga4_i
      use diaga5_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ns, ja1, ja2, k1, k2, ka, ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij, iskr, npeelgg
      real(double) :: s, ss, re
!-----------------------------------------------
      iat=1
      ij=jlist(ja1)
      s=dble(jjq1(3,ij))
      ij=jlist(ja2)
      ss=dble(jjq1(3,ij))
      ss=s*ss
      recc=one/dsqrt(ss)
      if(npeel == 1 .and. ns == -1)return
      if(ns == -1) then
         npeelgg = npeel
      else
         npeelgg = ns
      end if
      iat=0
      iskr=npeelgg-ja2
      if(iskr > 1) then
        call diaga3(ja2,npeelgg,2*ka,ire,iat,re)
        if(iat == 0)return
        recc=re*recc
        iat=0
      end if
      if(ja2 /= npeelgg) then
        call diaga5(npeelgg,ja2,2*ka,ire,iat,re)
        if(iat == 0)return
        recc=re*recc
        iat=0
      endif
      call diaga4(ja1,ja2,k1,k2,2*ka,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      if(ja1 == 1.and.ja2 == 2)return
      iat=0
      call diaga1(ja1,k1,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      iskr=ja2-ja1
      if(ja1 == 1)iskr=ja2-1-ja1
      if(iskr <= 1)return
      iat=0
      call diaga3(ja1,ja2,k1,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      return
      end subroutine recop2

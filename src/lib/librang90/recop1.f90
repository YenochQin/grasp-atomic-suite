!*******************************************************************
!                                                                  *
      subroutine recop1(ns,ja1,ka,ire,iat,recc)
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
      use m_c,             only: jlist, jjq1, npeel
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use diaga1_i
      use diaga3_i
      use diaga5_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ns, ja1, ka, ire
      integer, intent(out)      :: iat
      real(double), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ij1, iskr, npeelgg
      real(double) :: s, re
!-----------------------------------------------
      iat=1
      ij1=jlist(ja1)
      s=dble(jjq1(3,ij1))
      recc=one/dsqrt(s)
      if(npeel == 1 .and. ns == -1)return
      if(ns == -1) then
         npeelgg = npeel
      else
         npeelgg = ns
      end if
      iat=0
      if(ire /= 0) then
        if(ka == 0) then
          iat=1
          return
        end if
      end if
      iat=1
      if(npeelgg == 1) return
      iat=0
      if(npeelgg /= 2) then
        call diaga5(npeelgg,ja1,2*ka,ire,iat,re)
        recc=re*recc
        if(iat == 0) return
        if(ja1 == npeelgg) return
        iat=0
      end if
      call diaga1(ja1,2*ka,ire,iat,re)
      if(iat == 0)return
      recc=re*recc
      if(npeelgg == 2) return
      iskr=npeelgg-ja1
      if(ja1 == 1)iskr=npeelgg-1-ja1
      if(iskr <= 1)return
      iat=0
      call diaga3(ja1,npeelgg,2*ka,ire,iat,re)
      recc=re*recc
      return
      end subroutine recop1

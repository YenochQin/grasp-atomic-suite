!*******************************************************************
!                                                                  *
      subroutine el32(jjja,jjjb,ja,jb,jja,jjb,jjc,jjd,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 07  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2221, 2212 ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1222,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,RECO2,SIXJ,SPEAK                     *
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
      use cons_c,          only: zero, half, eps
      use m_c,             only: nq1, jlist, npeel
      use orb_c,           only: nak
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reco_i
      use reco2_i
      use perko2_i
      use coulom_i
      use snrc_i
      use ixjtik_i
      use itrexg_i
      use gg1112_i
      use sixj_i
      use speak_i
      use cxk_i
      use talk_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jjja,jjjb,ja,jb,jja,jjb,jjc,jjd,icolbrei
!      DIMENSION J(2)
!      DIMENSION S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iia,iib,iic,iid,iatikk,ip1,ig1,ii,i2,i3,ibrd,ibre,    &
                 ia,ib,iat,ifaz,ifazfrcs,ikk,inn,jb1,jaa,jbb,j12,kra,  &
                 l1,l2,n,nn,nd1,nd2,ne1,ne2
      integer, dimension(2) :: j
      integer, dimension(4) :: is,kaps,ks
      real(double)          :: aa,ab,a1,bb,qm1,qm2,qm3,qm4,si,recc
      real(double), dimension(12) :: s
!-----------------------------------------------
      if(npeel <= 1)return
      iia=jlist(jja)
      iib=jlist(jjb)
      iic=jlist(jjc)
      iid=jlist(jjd)
      if(ja > jb) then
        jaa=jb
        jbb=ja
      else
        jaa=ja
        jbb=jb
      end if
      call reco(jaa,jbb,jbb,jbb,1,iat)
      if(iat == 0)return
      ia=jlist(ja)
      ib=jlist(jb)
      qm1=half
      qm2=half
      qm3=-half
      qm4=-half
      call perko2(ja,jb,ja,ja,2)
      j(1)=id1(3)
      j(2)=id2(3)
      l1=(j(1)+1)/2
      l2=(j(2)+1)/2
      call reco2(jaa,jbb,j(1),0,iat,recc)
      if(iat == 0)return
      ip1=itrexg(j(2),j(2),j(1),j(2),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
      if (icolbrei == 2) then
        is(1)=iia
        is(2)=iib
        is(3)=iic
        is(4)=iid
        kaps(1)=2*nak(is(1))
        kaps(2)=2*nak(is(2))
        kaps(3)=2*nak(is(3))
        kaps(4)=2*nak(is(4))
        ks(1)=iabs(kaps(1))
        ks(2)=iabs(kaps(2))
        ks(3)=iabs(kaps(3))
        ks(4)=iabs(kaps(4))
        call snrc(is,kaps,ks,nd1,nd2,ne1,ne2,ibrd,ibre)
        if(ibrd <= 0)return
      end if
      call reco2(jaa,jbb,j(1),1,iat,recc)
! * * *                      * * *                      * * *
!     CASES 2221   + + - -        TRANSFORM TO  1222   - + + -
!           2212                                1222
!
      do i2=ip1,ig1,2
        kra=(i2-1)/2
        if (icolbrei == 1) then
          call coulom(l2,l2,l2,l1,id2(5),id2(5),id2(5),id1(5),kra,a1)
          if(dabs(a1) < eps) cycle
          a1=-a1
        end if
!
        ab=zero
        do i3=ip1,ig1,2
          j12=(i3-1)/2
          ifaz=j(1)-j12+1
          if((ifaz/2)*2 /= ifaz) cycle
          if(ixjtik(j(2),j(2),kra*2,j(1),j(2),j12*2) == 0) cycle
          call gg1112(ik2,ik1,bk2,bk1,id2,id1,bd2,     &
                      bd1,j12,qm1,qm2,qm3,qm4,aa)
          if(dabs(aa) < eps) cycle
          call sixj(j(2),j(2),kra*2,j(1),j(2),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ifaz=2*j(2)+kra*2+j12*2
          if((ifaz/4)*4 /= ifaz)aa=-aa
          ab=ab+aa
        end do
        ab=ab*recc
        if(dabs(ab) < eps) cycle
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        ifazfrcs=1
        ifaz=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)
        if((ifaz/4)*4 /= ifaz)ifazfrcs=-ifazfrcs
!
        nn=0
        jb1=jbb-1
        do ii=jaa,jb1
          inn=jlist(ii)
          nn=nq1(inn)+nn
        end do
        if((nn/2)*2 == nn)ab=-ab
        if (icolbrei == 1) then
          bb=ab*a1*dble(ifazfrcs)
          call speak(jjja,jjjb,iia,iib,iic,iid,kra,bb)
        else if (icolbrei == 2) then
          n=(kra-nd1)/2+1
          if(((kra-nd1)/2)*2 == (kra-nd1)) then
            call cxk(s,is,kaps,kra,kra,ibrd,1)
            if(dabs(s(1)) > eps) then
              bb=-s(1)*ab
              if(dabs(bb) > eps)                              &
              call talk(jjja,jjjb,kra,is(1),is(3),is(2),is(4),3,bb)
            end if
          end if
        end if
      end do
      return
      end subroutine el32

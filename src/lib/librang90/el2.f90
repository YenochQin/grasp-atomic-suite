!*******************************************************************
!                                                                  *
      subroutine el2(jja,jjb,ja,jb,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
!                                              N'2 = N2 + 2        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,       *
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
      use m_c,             only: jlist, npeel
      use orb_c,           only: nak
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reco_i
      use reco2_i
      use itrexg_i
      use itrig_i
      use ixjtik_i
      use snrc_i
      use coulom_i
      use gg1122_i
      use sixj_i
      use cxk_i
      use talk_i
      use perko2_i, only: perko2
      use speak_i, only: speak
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,icolbrei
!      DIMENSION J(2)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iat,ia,ib,ibrd,ibre,ip1,ip2,ig1,ig2,ikk,ii,i2,i3, &
                 ifaz,ifaz1,ifazfrcs,j12,jaa,jbb, &
                 kra,l1,l2,n,nd1,nd2,ne1,ne2,nup1,nu,mu
      integer, dimension(2) :: j
      integer, dimension(4) :: is,kaps,ks
      real(double)          :: aa,a1,ab,bb,qm1,qm2,qm3,qm4,recc,si
      real(double), dimension(12)    :: s
      real(double), dimension(12,20) :: cond
!-----------------------------------------------
      if(npeel <= 1)return
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
      ip1=itrexg(j(1),j(1),j(2),j(2),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
! * * *                      * * *                      * * *
!     THE CASE 1122   + + - -
!
      ip2=itrexg(j(1),j(2),j(1),j(2),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      if (icolbrei == 2) then
        is(1)=ia
        is(2)=ia
        is(3)=ib
        is(4)=ib
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
        do ii=1,20
          cond(1,ii) =zero
          cond(2,ii) =zero
          cond(3,ii) =zero
          cond(4,ii) =zero
          cond(5,ii) =zero
          cond(6,ii) =zero
          cond(7,ii) =zero
          cond(8,ii) =zero
          cond(9,ii) =zero
          cond(10,ii)=zero
          cond(11,ii)=zero
          cond(12,ii)=zero
        end do
      end if
      do i2=ip2,ig2,2
        kra=(i2-1)/2
        if (icolbrei == 1) then
          call coulom(l1,l1,l2,l2,id1(5),id1(5),id2(5),id2(5),kra,a1)
          if(dabs(a1) < eps) cycle
          a1=-half*a1
        end if
        ab=zero
          do i3=ip1,ig1,2
            j12=(i3-1)/2
            call reco2(jaa,jbb,j12*2,0,iat,recc)
            if(iat /= 0) then
              if(ixjtik(j(1),j(2),kra*2,j(2),j(1),j12*2) /= 0) then
                call gg1122(j12,j12,qm1,qm2,qm3,qm4,aa)
                if(dabs(aa) > eps) then
                  call reco2(jaa,jbb,j12*2,1,iat,recc)
                  aa=aa*recc
                  call sixj(j(1),j(2),kra*2,j(2),j(1),j12*2,0,si)
                  aa=aa*si*dsqrt(dble(i3))
                  ifaz=ik1(3)+ik2(3)+kra*2+j12*2
                  if((ifaz/4)*4 /= ifaz)aa=-aa
                  ab=ab+aa
                end if
              end if
            end if
          end do
!
!       TRANSFORM FANO & RACAH PHASE CONVENTION
!       TO CONDON & SHORTLEY PHASE CONVENTION
!
        ifazfrcs=1
        ifaz1=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)
        if((ifaz1/4)*4 /= ifaz1)ifazfrcs=-ifazfrcs
!
        if (icolbrei == 1) then
          bb=a1*ab*dble(ifazfrcs)
          if(dabs(bb) > eps)call speak(jja,jjb,ia,ia,ib,ib,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu+1) /= 0) .and.   &
               (itrig(ks(2),ks(4),nu+nu+1) /= 0)) then
              n=(nu-nd1)/2+1
              if(nu > 0) then
                call cxk(s,is,kaps,nu,kra,1,1)
                do mu = 1,4
                  cond(mu,n)=cond(mu,n)-half*ab*s(mu)
                end do
              end if
            end if
          end if
          nu=kra+1
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu-1) /= 0) .and.  &
               (itrig(ks(2),ks(4),nu+nu-1) /= 0)) then
              n=(nu-nd1)/2+1
              if(n <= nd2) then
                call cxk(s,is,kaps,nu,kra,1,1)
                do mu = 1,4
                  cond(mu,n)=cond(mu,n)-half*ab*s(mu)
                end do
              end if
            end if
          end if
          nu=kra-1
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu+3) /= 0) .and.   &
               (itrig(ks(2),ks(4),nu+nu+3) /= 0)) then
              if(nu >= 0) then
                n=(nu-nd1)/2+1
                if(n < nd2) then
                  call cxk(s,is,kaps,nu,kra,1,1)
                  do mu = 1,12
                    cond(mu,n)=cond(mu,n)-half*ab*s(mu)
                  end do
                end if
              end if
            end if
          end if
        end if
      end do
      if (icolbrei .eq. 2) then
        do n = 1,nd2
          nu=nd1+2*(n-1)
          call talk(jja,jjb,nu,is(1),is(3),is(2),is(4),1,cond(1,n))
          call talk(jja,jjb,nu,is(3),is(1),is(4),is(2),1,cond(2,n))
          call talk(jja,jjb,nu,is(1),is(3),is(4),is(2),1,cond(3,n))
          call talk(jja,jjb,nu,is(3),is(1),is(2),is(4),1,cond(4,n))
          if(n == nd2) cycle
          nup1=nu+1
          call talk(jja,jjb,nup1,is(1),is(3),is(2),is(4),2,cond(5,n))
          call talk(jja,jjb,nup1,is(2),is(4),is(1),is(3),2,cond(6,n))
          call talk(jja,jjb,nup1,is(3),is(1),is(4),is(2),2,cond(7,n))
          call talk(jja,jjb,nup1,is(4),is(2),is(3),is(1),2,cond(8,n))
          call talk(jja,jjb,nup1,is(1),is(3),is(4),is(2),2,cond(9,n))
          call talk(jja,jjb,nup1,is(4),is(2),is(1),is(3),2,cond(10,n))
          call talk(jja,jjb,nup1,is(3),is(1),is(2),is(4),2,cond(11,n))
          call talk(jja,jjb,nup1,is(2),is(4),is(3),is(1),2,cond(12,n))
        end do
      end if
      return
      end subroutine el2

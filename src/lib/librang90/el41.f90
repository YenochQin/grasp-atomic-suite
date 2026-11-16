!*******************************************************************
!                                                                  *
      subroutine el41(jjja,jjjb,ja,jb,jc,irez,jja,jjb,jjc,jjd,    &
                                                          icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 10  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
!                                              N'2 = N2 + 1        *
!                                              N'3 = N3 - 2,       *
!     WHEN IREZ = 1   . . . . . . . . . . . . . . . . . . .        *
!                                              N'1 = N1 - 1        *
!                                              N'2 = N2 - 1        *
!                                              N'3 = N3 + 2,       *
!     WHEN IREZ = 2   . . . . . . . . . . . . . . . . . . .        *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
!                        JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK         *
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
      use orb_c,           only: orbital_kappa
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use eile_i
      use reco_i
      use rec3_i
      use perko2_i
      use snrc_i
      use jfaze_i
      use ixjtik_i
      use itrexg_i
      use itrig_i
      use coulom_i
      use gg1233_i
      use sixj_i
      use speak_i
      use cxk_i
      use talk_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jjja,jjjb,ja,jb,jc,irez,jja,jjb, &
                             jjc,jjd,icolbrei
!      DIMENSION J(3)
!      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  ia,ib,ic,ii,iia,iib,iic,ibrd,ibre,ip1,ig1,ip2,ig2, &
                  iat,iid,ikk,ifaz,ifazp,ifazfrcs,inn,i2,i3,jb1,jaa, &
                  jbb,jcc,j12,kra,l1,l2,l3,nd1,nd2,ne1,ne2,n,nn,nu,  &
                  nup1,mu
      integer, dimension(3) :: j
      integer, dimension(4) :: is,kaps,ks
      real(double)          :: aa,ab,a1,bb,qm1,qm2,qm3,qm4,si,recc
      real(double), dimension(12) :: s
      real(double), dimension(12,20) :: cond
!-----------------------------------------------
      call eile(ja,jb,jc,jaa,jbb,jcc)
      if(npeel <= 1)return
      call reco(jaa,jcc,jbb,jbb,2,iat)
      if(iat == 0)return
      ia=jlist(ja)
      ib=jlist(jb)
      ic=jlist(jc)
      iia=jlist(jja)
      iib=jlist(jjb)
      iic=jlist(jjc)
      iid=jlist(jjd)
      if(irez == 1) then
        qm1=-half
        qm2=-half
        qm3=half
        qm4=half
      else
        qm1=half
        qm2=half
        qm3=-half
        qm4=-half
      end if
      call perko2(ja,jb,jc,ja,3)
      j(1)=id1(3)
      j(2)=id2(3)
      j(3)=id3(3)
      l1=(j(1)+1)/2
      l2=(j(2)+1)/2
      l3=(j(3)+1)/2
      if (icolbrei == 2) then
        is(1)=iia
        is(2)=iib
        is(3)=iic
        is(4)=iid
        kaps(1)=2*orbital_kappa(is(1))
        kaps(2)=2*orbital_kappa(is(2))
        kaps(3)=2*orbital_kappa(is(3))
        kaps(4)=2*orbital_kappa(is(4))
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
      ifazp=jfaze(jc,ja,jb,jc)
      ifazfrcs = 1
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      ifaz=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)+ &
      ik3(5)*ik3(4)-id3(5)*id3(4)
      if((ifaz/4)*4 /= ifaz)ifazfrcs=-ifazfrcs
!
      if(ja > jb) then
        jaa=jb
        jbb=ja
      else
        jaa=ja
        jbb=jb
      end if
      nn=0
      jb1=jbb-1
      do ii=jaa,jb1
        inn=jlist(ii)
        nn=nq1(inn)+nn
      end do
      if((nn/2)*2 == nn)ifazp=-ifazp
! * * *                      * * *                      * * *
!     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
!           3321                                1233
!                                                    (IREZ = 1)
!     OR
!     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
!           2133                                1233
!                                                    (IREZ = 2)
      ip1=itrexg(j(2),j(1),j(3),j(3),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
      ip2=itrexg(j(3),j(1),j(2),j(3),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      do i2=ip2,ig2,2
        kra=(i2-1)/2
!
        if (icolbrei == 1) then
          if(irez == 2) then
            call coulom(l1,l2,l3,l3,id1(5),id2(5),id3(5),id3(5),kra,a1)
          else
            call coulom(l3,l3,l1,l2,id3(5),id3(5),id1(5),id2(5),kra,a1)
          end if
          if(dabs(a1) < eps) cycle
        end if
!
        ab=zero
        do i3=ip1,ig1,2
          j12=(i3-1)/2
          ifaz=j(2)-j12+1
          if(irez == 2)ifaz=j(1)-j12+1
          if((ifaz/2)*2 /= ifaz) cycle
          call rec3(ja,jb,jc,j(1),j(2),j12*2,0,iat,aa)
          if(iat == 0) cycle
          if(ixjtik(j(3),j(1),kra*2,j(2),j(3),j12*2) == 0) cycle
          call gg1233(ik1,ik2,ik3,bk1,bk2,bk3,id1,id2,id3,bd1,    &
                    bd2,bd3,j12,qm1,qm2,qm3,qm4,aa)
          if(dabs(aa) < eps) cycle
          call rec3(ja,jb,jc,j(1),j(2),j12*2,1,iat,recc)
          aa=aa*recc
          call sixj(j(3),j(1),kra*2,j(2),j(3),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ifaz=j(3)+j(1)+2*j12+2*kra
          if(irez == 2)ifaz=j(2)+j(3)+2*j12+2*kra
          if((ifaz/4)*4 /= ifaz)aa=-aa
          ab=ab+aa
        end do
        if(dabs(ab) < eps) cycle
        ab=-ab*dble(ifazp)
        if (icolbrei == 1) then
          bb=a1*ab*dble(ifazfrcs)
          call speak(jjja,jjjb,iia,iib,iic,iid,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu+1) /= 0) .and.  &
               (itrig(ks(2),ks(4),nu+nu+1) /= 0)) then
              n=(nu-nd1)/2+1
              if(nu > 0) then
                call cxk(s,is,kaps,nu,kra,1,1)
                do mu = 1,4
                  cond(mu,n)=cond(mu,n)+ab*s(mu)
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
                  cond(mu,n)=cond(mu,n)+ab*s(mu)
                end do
              end if
            end if
          end if
          nu=kra-1
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu+3) /= 0) .and.  &
               (itrig(ks(2),ks(4),nu+nu+3) /= 0)) then
              if(nu >= 0) then
                n=(nu-nd1)/2+1
                if(n < nd2) then
                  call cxk(s,is,kaps,nu,kra,1,1)
                  do mu = 1,12
                    cond(mu,n)=cond(mu,n)+ab*s(mu)
                  end do
                end if
              end if
            end if
          end if
        end if
      end do
      if (icolbrei == 2) then
        do n = 1,nd2
          nu=nd1+2*(n-1)
          call talk(jjja,jjjb,nu,is(1),is(3),is(2),is(4),1,cond(1,n))
          call talk(jjja,jjjb,nu,is(3),is(1),is(4),is(2),1,cond(2,n))
          call talk(jjja,jjjb,nu,is(1),is(3),is(4),is(2),1,cond(3,n))
          call talk(jjja,jjjb,nu,is(3),is(1),is(2),is(4),1,cond(4,n))
          if(n == nd2) cycle
          nup1=nu+1
          call talk(jjja,jjjb,nup1,is(1),is(3),is(2),is(4),2,cond(5,n))
          call talk(jjja,jjjb,nup1,is(2),is(4),is(1),is(3),2,cond(6,n))
          call talk(jjja,jjjb,nup1,is(3),is(1),is(4),is(2),2,cond(7,n))
          call talk(jjja,jjjb,nup1,is(4),is(2),is(3),is(1),2,cond(8,n))
          call talk(jjja,jjjb,nup1,is(1),is(3),is(4),is(2),2,cond(9,n))
          call talk(jjja,jjjb,nup1,is(4),is(2),is(1),is(3),2,cond(10,n))
          call talk(jjja,jjjb,nup1,is(3),is(1),is(2),is(4),2,cond(11,n))
          call talk(jjja,jjjb,nup1,is(2),is(4),is(3),is(1),2,cond(12,n))
        end do
      end if
      return
   10 write(99,100)
  100 format(5x,'ERRO IN EL41  PMGG RAGG')
      stop
      end subroutine el41

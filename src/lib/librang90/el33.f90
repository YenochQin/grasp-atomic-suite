module el33_i
    implicit none (type, external)
    private
    public :: el33
contains
!*******************************************************************
!                                                                  *
      subroutine el33(jjja,jjjb,ja,jb,jc,irez,jja,jjb,jjc,jjd   &
                                                        ,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 08  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
!                      JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK           *
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
      use itrexg_i
      use itrig_i
      use ixjtik_i
      use jfaze_i
      use snrc_i
      use gg1233_i
      use coulom_i
      use speak_i
      use cxk_i
      use talk_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jjja,jjjb,ja,jb,jc,irez,jja,jjb,jjc,jjd &
                                                        ,icolbrei
!      DIMENSION PMGG(30),RAGG(30),J(3)
!      DIMENSION CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ia,ib,ic,iat,ii,iia,iib,iic,iid,ip1,ip2,ig1,ig2,   &
                 ikk,inn,i1,i2,i3,i4,ibrd,ibre,ifaz,ifazp,ifazfrcs, &
                 jaa,jbb,jcc,jb1,j12,kra,kra1,l1,l2,l3,nd1,nd2,ne1, &
                 ne2,n,nn,nu,nup1,mu
      integer :: j(3)
      integer :: is(4),kaps(4),ks(4)
      real(real64)          :: a1,aa,ab,bb,qm1,qm2,qm3,qm4,rag,recc,si
      real(real64) :: s(12)
      real(real64) :: pmgg(30),ragg(30)
      real(real64) :: cone(12, 20)
!-----------------------------------------------
      if(npeel <= 1)return
      call eile(ja,jb,jc,jaa,jbb,jcc)
      call reco(jaa,jcc,jbb,jbb,2,iat)
      if(iat == 0)return
      ia=jlist(ja)
      ib=jlist(jb)
      ic=jlist(jc)
      iia=jlist(jja)
      iib=jlist(jjb)
      iic=jlist(jjc)
      iid=jlist(jjd)
      qm1=half
      qm2=-half
      qm3=half
      qm4=-half
      call perko2(ja,jb,jc,ja,3)
      j(1)=id1(3)
      j(2)=id2(3)
      j(3)=id3(3)
      l1=(j(1)+1)/2
      l2=(j(2)+1)/2
      l3=(j(3)+1)/2
      ip1=itrexg(j(2),j(1),j(3),j(3),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
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
        if(ibrd <= 0 .and. ibre <= 0)return
        do ii=1,20
          cone(1,ii) =zero
          cone(2,ii) =zero
          cone(3,ii) =zero
          cone(4,ii) =zero
          cone(5,ii) =zero
          cone(6,ii) =zero
          cone(7,ii) =zero
          cone(8,ii) =zero
          cone(9,ii) =zero
          cone(10,ii)=zero
          cone(11,ii)=zero
          cone(12,ii)=zero
        end do
      end if
      do i4=ip1,ig1,2
        kra=(i4-1)/2
        kra1=kra+1
        if(kra1 > 30)go to 10
        ragg(kra1)=zero
        pmgg(kra1)=zero
        call rec3(jb,ja,jc,j(2),j(1),kra*2,0,iat,recc)
        if(iat == 0) cycle
        call gg1233(ik2,ik1,ik3,bk2,bk1,bk3,id2,id1,id3,bd2,  &
                    bd1,bd3,kra,qm1,qm2,qm3,qm4,rag)
        if(dabs(rag) < eps) cycle
        ragg(kra1)=rag
        call rec3(jb,ja,jc,j(2),j(1),kra*2,1,iat,recc)
        pmgg(kra1)=recc
      end do
      ifazp=jfaze(jb,ja,jc,jc)
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      ifazfrcs=1
      ifaz=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)
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
      if(irez == 2)go to 5
! * * *                      * * *                      * * *
!     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
!           3231                                2133
!
    6 continue
      do i1=ip1,ig1,2
        kra=(i1-1)/2
        kra1=kra+1
        if(kra1 > 30)go to 10
        if (icolbrei == 1) then
          call coulom(l2,l3,l1,l3,id2(5),id3(5),id1(5),id3(5),kra,a1)
          if(dabs(a1) < eps) cycle
        end if
        aa=pmgg(kra1)
        if(dabs(aa) < eps) cycle
        aa=aa*ragg(kra1)
        if(dabs(aa) < eps) cycle
        aa=aa/dsqrt(dble(i1))
        aa=aa*dble(ifazp)
        if (icolbrei == 1) then
          bb=a1*aa*dble(ifazfrcs)
          call speak(jjja,jjjb,iia,iib,iic,iid,kra,bb)
        else if (icolbrei == 2) then
          n=(kra-nd1)/2+1
          call cxk(s,is,kaps,kra,kra,2,1)
          if(dabs(s(1)) > eps) then
            bb=s(1)*aa
            if(dabs(bb) > eps)                    &
          call talk(jjja,jjjb,kra,is(1),is(3),is(2),is(4),3,bb)
          end if
        end if
      end do
      if(irez == 2) go to 7
! * * *                      * * *                      * * *
!     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
!           2331                                2133
!
    5 ip2=itrexg(j(3),j(1),j(2),j(3),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      do i2=ip2,ig2,2
        kra=(i2-1)/2
        if(kra > 30)go to 10
        if (icolbrei == 1) then
          call coulom(l3,l2,l1,l3,id3(5),id2(5),id1(5),id3(5),kra,a1)
          if(dabs(a1) < eps) cycle
        end if
        ab=zero
        do i3=ip1,ig1,2
          j12=(i3-1)/2
          kra1=j12+1
          if(kra1 > 30)go to 10
          aa=pmgg(kra1)
          if(dabs(aa) < eps) cycle
          aa=aa*ragg(kra1)
          if(dabs(aa) < eps) cycle
          if(ixjtik(j(1),j(3),kra*2,j(3),j(2),j12*2) == 0) cycle
          call sixj(j(1),j(3),kra*2,j(3),j(2),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ab=ab+aa
        end do
        if(dabs(ab) < eps) cycle
        ab=ab*dble(ifazp)
        if (icolbrei == 1) then
          bb=a1*ab*dble(ifazfrcs)
          call speak(jjja,jjjb,iia,iib,iid,iic,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu+1) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu+1) /= 0)) then
              if(nu > 0) then
                n=(nu-ne1)/2+1
                call cxk(s,is,kaps,nu,kra,1,2)
                do mu = 1,4
                  cone(mu,n)=cone(mu,n)+ab*s(mu)
                end do
              end if
            end if
          end if
          nu=kra+1
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu-1) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu-1) /= 0)) then
              if(nu >= 0) then
                n=(nu-ne1)/2+1
                if(n <= ne2) then
                  call cxk(s,is,kaps,nu,kra,1,2)
                  do mu = 1,4
                    cone(mu,n)=cone(mu,n)+ab*s(mu)
                  end do
                end if
              end if
            end if
          end if
          nu=kra-1
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu+3) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu+3) /= 0)) then
              if(nu >=  0) then
                n=(nu-ne1)/2+1
                if(n < ne2) then
                  call cxk(s,is,kaps,nu,kra,1,2)
                  do mu = 1,12
                    cone(mu,n)=cone(mu,n)+ab*s(mu)
                  end do
                end if
              end if
            end if
          end if
        end if
      end do
      if (icolbrei == 2) then
        do n = 1,ne2
          nu=ne1+2*(n-1)
          call talk(jjja,jjjb,nu,is(1),is(4),is(2),is(3),1,cone(1,n))
          call talk(jjja,jjjb,nu,is(4),is(1),is(3),is(2),1,cone(2,n))
          call talk(jjja,jjjb,nu,is(1),is(4),is(3),is(2),1,cone(3,n))
          call talk(jjja,jjjb,nu,is(4),is(1),is(2),is(3),1,cone(4,n))
          if(n == ne2) cycle
          nup1=nu+1
          call talk(jjja,jjjb,nup1,is(1),is(4),is(2),is(3),2,cone(5,n))
          call talk(jjja,jjjb,nup1,is(2),is(3),is(1),is(4),2,cone(6,n))
          call talk(jjja,jjjb,nup1,is(4),is(1),is(3),is(2),2,cone(7,n))
          call talk(jjja,jjjb,nup1,is(3),is(2),is(4),is(1),2,cone(8,n))
          call talk(jjja,jjjb,nup1,is(1),is(4),is(3),is(2),2,cone(9,n))
          call talk(jjja,jjjb,nup1,is(3),is(2),is(1),is(4),2,cone(10,n))
          call talk(jjja,jjjb,nup1,is(4),is(1),is(2),is(3),2,cone(11,n))
          call talk(jjja,jjjb,nup1,is(2),is(3),is(4),is(1),2,cone(12,n))
        end do
      end if
      if(irez == 2)go to 6
    7 continue
      return
   10 write(99,100)
  100 format(5x,'ERRO IN EL33  PMGG RAGG')
      stop
      end subroutine el33
end module el33_i

module el52_i
    implicit none (type, external)
    private
    public :: el52
contains
!*******************************************************************
!                                                                  *
      subroutine el52(jja,jjb,ja,jb,jc,jd,irez,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 13  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 1324, 3142, 1342, 3124    *
!                                                   ( IREZ = 1),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 - 1   *
!                                                   N'2 = N2 + 1   *
!                                                   N'3 = N3 - 1   *
!                                                   N'4 = N4 + 1   *
!     AND    2413, 4231, 2431, 4213                 ( IREZ = 2),   *
!                                                   ( + + - - ),   *
!     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
!     CONFIGURATIONS:                               N'1 = N1 + 1   *
!                                                   N'2 = N2 - 1   *
!                                                   N'3 = N3 + 1   *
!                                                   N'4 = N4 - 1   *
!                                                                  *
!     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
!                        RECO,REC4,SIXJ,SPEAK                      *
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
      use reco_i
      use reco4_i
      use perko2_i
      use snrc_i
      use gg1234_i
      use itrexg_i
      use ixjtik_i
      use coulom_i
      use sixj_i
      use speak_i
      use itrig_i
      use cxk_i
      use talk_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,jc,jd,irez,icolbrei
!      DIMENSION PMGG(30),J(4)
!     DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ia,ib,ic,id,ibrd,ibre,ii,ip1,ip2,ig1,ig2,ikk,i2,i3,  &
                 i4,ifaz,ifazp,ifazfrcs,inn,iat,kra,kra1,l1,l2,l3,l4, &
                 j12,jb1,jd1,nd1,nd2,ne1,ne2,n,nn,nu,nup1,mu
      integer :: j(4)
      integer :: is(4),kaps(4),ks(4)
      real(real64)          :: aa,ab,a1,bb,qm1,qm2,qm3,qm4,rag,recc,si
      real(real64) :: s(12)
      real(real64) :: pmgg(30)
      real(real64) :: cond(12, 20),cone(12, 20)
!-----------------------------------------------
      if(npeel <= 3)return
      call reco(ja,jd,jc,jb,3,iat)
      if(iat == 0)return
      ia=jlist(ja)
      ib=jlist(jb)
      ic=jlist(jc)
      id=jlist(jd)
      if(irez == 2) then
        qm1=-half
        qm2=half
        qm3=-half
        qm4=half
      else
        qm1=half
        qm2=-half
        qm3=half
        qm4=-half
      end if
      call perko2(ja,jb,jc,jd,4)
      j(1)=id1(3)
      j(2)=id2(3)
      j(3)=id3(3)
      j(4)=id4(3)
      l1=(j(1)+1)/2
      l2=(j(2)+1)/2
      l3=(j(3)+1)/2
      l4=(j(4)+1)/2
      if (icolbrei == 2) then
        if(irez == 1) then
          is(1)=ia
          is(2)=ic
          is(3)=ib
          is(4)=id
        else if(irez == 2) then
          is(1)=ib
          is(2)=id
          is(3)=ia
          is(4)=ic
        end if
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
      call gg1234(ik1,ik2,ik3,ik4,bk1,bk2,bk3,bk4,id1,id2,  &
      id3,id4,bd1,bd2,bd3,bd4,qm1,qm2,qm3,qm4,rag)
      if(dabs(rag) < eps) return
      ip1=itrexg(j(1),j(2),j(3),j(4),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
      do i4=ip1,ig1,2
        kra=(i4-1)/2
        kra1=kra+1
        if(kra1 > 30)go to 10
        pmgg(kra1)=zero
        call reco4(ja,jb,jc,jd,j(1),j(2),j(3),j(4),kra*2,0,iat,recc)
        if(iat == 0) cycle
        call reco4(ja,jb,jc,jd,j(1),j(2),j(3),j(4),kra*2,1,iat,recc)
        pmgg(kra1)=recc
      end do
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
      ifazfrcs=1
      ifaz=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)  &
      +ik3(5)*ik3(4)-id3(5)*id3(4)+ik4(5)*ik4(4)-id4(5)*id4(4)
      if((ifaz/4)*4 /= ifaz)ifazfrcs=-ifazfrcs
!
      nn=0
      jb1=jb-1
      ifazp=1
      do ii=ja,jb1
        inn=jlist(ii)
        nn=nq1(inn)+nn
      end do
      if((nn/2)*2 == nn)ifazp=-ifazp
      nn=0
      jd1=jd-1
      do ii=jc,jd1
        inn=jlist(ii)
        nn=nq1(inn)+nn
      end do
      if((nn/2)*2 == nn)ifazp=-ifazp
! * * *                      * * *                      * * *
!     CASES 1324   + + - -        TRANSFORM TO  1234   + - + -
!           1342                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2413   + + - -        TRANSFORM TO  1234   - + - +
!           4231                                1234
!                                                    (IREZ = 2)
      do i3=ip1,ig1,2
        kra=(i3-1)/2
        if (icolbrei == 1) then
          if(irez == 2) then
            call coulom(l2,l4,l1,l3,id2(5),id4(5),id1(5),id3(5),kra,a1)
          else
            call coulom(l1,l3,l2,l4,id1(5),id3(5),id2(5),id4(5),kra,a1)
          end if
          if(dabs(a1) < eps) cycle
        end if
        kra1=kra+1
        if(kra1 > 30)go to 10
        aa=pmgg(kra1)
        if(dabs(aa) < eps) cycle
        aa=aa*rag
        if(dabs(aa) < eps) cycle
        aa=aa/dsqrt(dble(i3))
        ab=aa*dble(ifazp)
        if(irez == 2) then
          ifaz=j(1)+j(2)+j(4)+j(3)+4*kra
          if((ifaz/4)*4 /= ifaz)ab=-ab
        end if
        if (icolbrei == 1) then
          bb=a1*ab*dble(ifazfrcs)
          if(irez == 1)call speak(jja,jjb,ia,ic,ib,id,kra,bb)
          if(irez == 2)call speak(jja,jjb,ib,id,ia,ic,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-nd1)/2)*2 == (nu-nd1)) then
            if((itrig(ks(1),ks(3),nu+nu+1) /= 0) .and.  &
               (itrig(ks(2),ks(4),nu+nu+1) /= 0)) then
              n=(nu-nd1)/2+1
              if(nu >  0) then
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
            if((itrig(ks(1),ks(3),nu+nu+3) /= 0) .and.   &
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
! * * *                      * * *                      * * *
!     CASES 1342   + + - -        TRANSFORM TO  1234   + - + -
!           3124                                1234
!                                                    (IREZ = 1)
!     OR
!     CASES 2431   + + - -        TRANSFORM TO  1234   - + - +
!           4213                                1234
!                                                    (IREZ = 2)
      ip2=itrexg(j(1),j(4),j(2),j(3),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      do i2=ip2,ig2,2
        kra=(i2-1)/2
!
        if (icolbrei == 1) then
          if(irez == 2) then
            call coulom(l2,l4,l3,l1,id2(5),id4(5),id3(5),id1(5),kra,a1)
          else
            call coulom(l1,l3,l4,l2,id1(5),id3(5),id4(5),id2(5),kra,a1)
          end if
          if(dabs(a1) < eps) cycle
        end if
!
        ab=zero
        do i3=ip1,ig1,2
          j12=(i3-1)/2
          ifaz=j(1)-j(3)+2*j12
          if(irez == 2)ifaz=j(4)-j(2)-2*j12
          if((ifaz/2)*2 /= ifaz) cycle
          kra1=j12+1
          if(kra1 > 30)go to 10
          aa=pmgg(kra1)
          if(dabs(aa) < eps) cycle
          aa=aa*rag
          if(dabs(aa) < eps) cycle
          if(ixjtik(j(1),j(4),kra*2,j(3),j(2),j12*2) == 0) cycle
          call sixj(j(1),j(4),kra*2,j(3),j(2),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ifaz=2*j(3)-4*kra-4*j12
          if(irez == 2)ifaz=j(1)+j(2)+j(3)-j(4)-4*kra
          if((ifaz/4)*4 /= ifaz)aa=-aa
          ab=ab+aa
        end do
        if(dabs(ab) < eps) cycle
        ab=-ab*dble(ifazp)
        if (icolbrei == 1) then
          bb=a1*ab*dble(ifazfrcs)
          if(irez == 1)call speak(jja,jjb,ia,ic,id,ib,kra,bb)
          if(irez == 2)call speak(jja,jjb,ib,id,ic,ia,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu+1) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu+1) /= 0)) then
              n=(nu-ne1)/2+1
              if(nu > 0) then
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
              n=(nu-ne1)/2+1
              if(n <= ne2) then
                call cxk(s,is,kaps,nu,kra,1,2)
                do mu = 1,4
                  cone(mu,n)=cone(mu,n)+ab*s(mu)
                end do
              end if
            end if
          end if
          nu=kra-1
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu+3) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu+3) /= 0)) then
              if(nu >= 0) then
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
          call talk(jja,jjb,nu,is(1),is(4),is(2),is(3),1,cone(1,n))
          call talk(jja,jjb,nu,is(4),is(1),is(3),is(2),1,cone(2,n))
          call talk(jja,jjb,nu,is(1),is(4),is(3),is(2),1,cone(3,n))
          call talk(jja,jjb,nu,is(4),is(1),is(2),is(3),1,cone(4,n))
          if(n == ne2) cycle
          nup1=nu+1
          call talk(jja,jjb,nup1,is(1),is(4),is(2),is(3),2,cone(5,n))
          call talk(jja,jjb,nup1,is(2),is(3),is(1),is(4),2,cone(6,n))
          call talk(jja,jjb,nup1,is(4),is(1),is(3),is(2),2,cone(7,n))
          call talk(jja,jjb,nup1,is(3),is(2),is(4),is(1),2,cone(8,n))
          call talk(jja,jjb,nup1,is(1),is(4),is(3),is(2),2,cone(9,n))
          call talk(jja,jjb,nup1,is(3),is(2),is(1),is(4),2,cone(10,n))
          call talk(jja,jjb,nup1,is(4),is(1),is(2),is(3),2,cone(11,n))
          call talk(jja,jjb,nup1,is(2),is(3),is(4),is(1),2,cone(12,n))
        end do
      end if
      return
   10 write(99,100)
  100 format(5x,'ERRO IN EL52  PMGG RAGG')
      stop
      end subroutine el52
end module el52_i

module el1_i
    implicit none (type, external)
    private
    public :: el1
contains
!*******************************************************************
!                                                                  *
      subroutine el1(jja,jjb,ja,jb,iire,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!      SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,      *
!                         RECO,RECO2,SIXJ,SPEAK,WW1                *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2020  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, half, eps
      use m_c,             only: jlist, npeel
      use orb_c,           only: orbital_kappa
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reco_i, only: reco
      use reco2_i, only: reco2
      use perko2_i, only: perko2
      use itrig_i, only: itrig
      use itrexg_i, only: itrexg
      use ixjtik_i, only: ixjtik
      use snrc_i, only: snrc
      use speak_i, only: speak
      use coulom_i, only: coulom
      use ww1_i, only: ww1
      use sixj_i, only: sixj
      use cxk_i, only: cxk
      use talk_i, only: talk
      use gg1122_i, only: gg1122
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,iire,icolbrei
!      DIMENSION CONE(7,20),S(12),IS(4),KAPS(4),KS(4)
!      DIMENSION PMGG(30),RAGG(30),J(2)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii,ia,ib,iat,ip1,ip2,ip3,ig1,ig2,ig3,ikk,i1,i2,i3,i4,&
                 ifaz,j12,ibrd,ibre,kra,kra1,l1,l2,mu,n,nu,nd1,nd2,   &
                 ne1,ne2,nup1
      integer :: j(2)
      integer :: is(4),kaps(4),ks(4)
      real(real64)          :: qm1,qm2,qm3,qm4,aa,ab,a1,bb,si,recc,rag
      real(real64)   :: s(12)
      real(real64)   :: pmgg(30),ragg(30)
      real(real64) :: cone(7, 20)
!-----------------------------------------------
      if(ja /= jb)goto 9
!
!     THE CASE 1111   + + - -
!
      if(iire /= 0) then
        call reco(ja,ja,ja,ja,0,iat)
        if(iat == 0)return
      end if
      call perko2(ja,ja,ja,ja,1)
      qm1=half
      qm2=half
      qm3=-half
      qm4=-half
      ia=jlist(ja)
      j(1)=id1(3)
      ip2=itrexg(j(1),j(1),j(1),j(1),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      l1=(j(1)+1)/2
      ip1=ip2
      ig1=ig2
      if (icolbrei == 2) then
        is(1)=ia
        is(2)=ia
        is(3)=ia
        is(4)=ia
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
      end if
      do i2=ip1,ig1,2
        kra=(i2-1)/2
        if (icolbrei == 1) then
          call coulom(l1,l1,l1,l1,id1(5),id1(5),id1(5),id1(5),kra,a1)
          if(dabs(a1) < eps)cycle
          a1=-a1*half
        end if
        ab=zero
        do i3=ip2,ig2,2
          j12=(i3-1)/2
          if(ixjtik(j(1),j(1),kra*2,j(1),j(1),j12*2) == 0)cycle
          call ww1(ik1,bk1,id1,bd1,j12,qm1,qm2,qm3,qm4,aa)
          if(dabs(aa) < eps)cycle
          call sixj(j(1),j(1),kra*2,j(1),j(1),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ifaz=ik1(3)+j12+kra
          if((ifaz/2)*2 /= ifaz)aa=-aa
          ab=ab+aa
        end do
!
!     RECOUPLING COEFFICIENTS
!
        if (icolbrei == 1) then
          bb=ab*a1
          bb=bb/dsqrt(dble(ik1(6)+1))
          if(dabs(bb) > eps)call speak(jja,jjb,ia,ia,ia,ia,kra,bb)
        else if (icolbrei == 2) then
          n=(kra-nd1)/2+1
          if(((kra-nd1)/2)*2 == (kra-nd1)) then
            call cxk(s,is,kaps,kra,kra,3,1)
            if(dabs(s(1)) > eps) then
              bb =-half*s(1)*ab/dsqrt(dble(ik1(6)+1))
              if(dabs(bb) > eps)call talk(jja,jjb,kra,ia,ia,ia,ia,4,bb)
            end if
          end if
        end if
      end do
      return
!  ............................................................
    9 if(npeel <= 1)return
      if(iire /= 0) then
        call reco(ja,jb,jb,jb,1,iat)
        if(iat == 0)return
      end if
      ia=jlist(ja)
      ib=jlist(jb)
      qm1=half
      qm2=-half
      qm3=half
      qm4=-half
      call perko2(ja,jb,ja,ja,2)
      j(1)=id1(3)
      j(2)=id2(3)
      l1=(j(1)+1)/2
      l2=(j(2)+1)/2
      ip1=itrexg(j(1),j(1),j(2),j(2),ikk)+1
      if(ikk <= 0)return
      ig1=ip1+ikk-1
      ip3=ip1
      ig3=ig1
      do i4=ip1,ig1,2
        kra=(i4-1)/2
        kra1=kra+1
        if(kra1 > 30)goto 10
        ragg(kra1)=zero
        pmgg(kra1)=zero
        call reco2(ja,jb,kra*2,0,iat,recc)
        if(iat == 0) cycle
        call gg1122(kra,kra,qm1,qm2,qm3,qm4,rag)
        if(dabs(rag) < eps) cycle
        ragg(kra1)=rag
        call reco2(ja,jb,kra*2,1,iat,recc)
        pmgg(kra1)=recc
      end do
! * * *                      * * *                      * * *
!     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
!           2121                                1122
!
      if (icolbrei == 2) then
        is(1)=ia
        is(2)=ib
        is(3)=ia
        is(4)=ib
        kaps(1)=2*orbital_kappa(is(1))
        kaps(2)=2*orbital_kappa(is(2))
        kaps(3)=2*orbital_kappa(is(3))
        kaps(4)=2*orbital_kappa(is(4))
        ks(1)=iabs(kaps(1))
        ks(2)=iabs(kaps(2))
        ks(3)=iabs(kaps(3))
        ks(4)=iabs(kaps(4))
        call snrc(is,kaps,ks,nd1,nd2,ne1,ne2,ibrd,ibre)
        do ii=1,20
           cone(1,ii)=zero
           cone(2,ii)=zero
           cone(3,ii)=zero
           cone(4,ii)=zero
           cone(5,ii)=zero
           cone(6,ii)=zero
           cone(7,ii)=zero
        end do
        if(ibrd == 0 .and. ibre == 0)return
      end if
      do i1=ip1,ig1,2
        kra=(i1-1)/2
        kra1=kra+1
        if(kra1 > 30)goto 10
        if (icolbrei == 1) then
          call coulom(l1,l2,l1,l2,id1(5),id2(5),id1(5),id2(5),kra,aa)
          if(dabs(aa) < eps) cycle
          aa=aa*pmgg(kra1)
          if(dabs(aa) < eps) cycle
          aa=aa*ragg(kra1)
          if(dabs(aa) < eps) cycle
          aa=aa/dsqrt(dble(i1))
          if(dabs(aa) > eps) call speak(jja,jjb,ia,ib,ia,ib,kra,aa)
        else if (icolbrei == 2) then
          n=(kra-nd1)/2+1
          if(((kra-nd1)/2)*2 == (kra-nd1)) then
            call cxk(s,is,kaps,kra,kra,3,1)
            if(dabs(s(1)) > eps) then
              bb=s(1)*pmgg(kra1)*ragg(kra1)/dsqrt(dble(i1))
              if(dabs(bb) > eps)call talk(jja,jjb,kra,ia,ia,ib,ib,4,bb)
            end if
          end if
        end if
      end do
! * * *                      * * *                      * * *
!     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
!           2112                                1122
!
      ip2=itrexg(j(1),j(2),j(1),j(2),ikk)+1
      if(ikk <= 0) return
      ig2=ip2+ikk-1
      do i2=ip2,ig2,2
        kra=(i2-1)/2
        if(kra > 30)goto 10
        if (icolbrei == 1) then
          call coulom(l1,l2,l2,l1,id1(5),id2(5),id2(5),id1(5),kra,a1)
          if(dabs(a1) < eps) cycle
        end if
        ab=zero
        do i3=ip3,ig3,2
          j12=(i3-1)/2
          kra1=j12+1
          if(kra1 > 30)goto 10
          aa=pmgg(kra1)
          if(dabs(aa) < eps) cycle
          aa=aa*ragg(kra1)
          if(dabs(aa) < eps) cycle
          if(ixjtik(j(1),j(2),kra*2,j(2),j(1),j12*2) == 0)cycle
          call sixj(j(1),j(2),kra*2,j(2),j(1),j12*2,0,si)
          aa=aa*si*dsqrt(dble(i3))
          ab=ab+aa
        end do
        if (icolbrei == 1) then
          bb=a1*ab
          if(dabs(bb) > eps)call speak(jja,jjb,ia,ib,ib,ia,kra,bb)
        else if (icolbrei == 2) then
          nu=kra
          if(((nu-ne1)/2)*2 == (nu-ne1)) then
            if((itrig(ks(1),ks(4),nu+nu+1) /= 0) .and.  &
               (itrig(ks(2),ks(3),nu+nu+1) /= 0)) then
              if(nu > 0) then
                n=(nu-ne1)/2+1
                call cxk(s,is,kaps,nu,kra,4,2)
                do mu = 1,3
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
                  call cxk(s,is,kaps,nu,kra,4,2)
                  do mu = 1,3
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
              if(nu >= 0) then
                n=(nu-ne1)/2+1
                if(n < ne2) then
                  call cxk(s,is,kaps,nu,kra,4,2)
                  do mu = 1,7
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
           call talk(jja,jjb,nu,ib,ia,ib,ia,5,cone(1,n))
           call talk(jja,jjb,nu,ia,ib,ib,ia,5,cone(2,n))
           call talk(jja,jjb,nu,ia,ib,ia,ib,5,cone(3,n))
           if(n == ne2) cycle
           nup1=nu+1
           call talk(jja,jjb,nup1,ia,ib,ia,ib,6,cone(4,n))
           call talk(jja,jjb,nup1,ib,ia,ib,ia,6,cone(5,n))
           call talk(jja,jjb,nup1,ia,ib,ib,ia,6,cone(6,n))
           call talk(jja,jjb,nup1,ib,ia,ia,ib,6,cone(7,n))
        end do
      end if
      return
   10 write(99,100)
  100 format(5x,"ERRO IN EL1  PMGG RAGG")
      stop
      end subroutine el1
end module el1_i

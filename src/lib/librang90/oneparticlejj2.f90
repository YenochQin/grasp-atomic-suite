      subroutine oneparticlejj2(ns,ka,ja,jb,coeff)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1 +- 1   *
!                                                  N'2 = N2 -+ 1   *
!                                                                  *
!      SUBROUTINE CALLED:                                          *
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
      use cons_c,          only: zero, tenth, half, eps
      use m_c,             only: nq1, jlist
      use orb_c,           only: nak
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use recop00_i
      use recop2_i
      use c0t5s_i
      use rmeajj_i
      use wj1_i
      use perko2_i, only: perko2
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ns,ka,ja,jb
      real(double), intent(out) :: coeff
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: i,ia,ib,iat,ij,ifaz,iqmm1,iqmm2,jibks1,jib, &
                      ks1,ks2
      real(double) :: rec,a1,a2,a3,s1,s2,qm1,qm2
!-----------------------------------------------
!
!     THE CASE 12   + -
!
      coeff=zero
      ia=min0(ja,jb)
      ib=max0(ja,jb)
      call recop00(ns,ia,ib,ka,iat)
      if(iat == 0)return
      ij=jlist(ia)
      ks1=(iabs(nak(ij))*2)-1
      ij=jlist(ib)
      ks2=(iabs(nak(ij))*2)-1
      call recop2(ns,ia,ib,ks1,ks2,ka,0,iat,rec)
      if(iat == 0)return
      call perko2(ja,jb,ja,ja,2)
      qm1=half
      qm2=-half
      iqmm1=qm1+qm1+tenth*qm1
      if(ik1(4) /= (id1(4)+iqmm1)) return
      iqmm2=qm2+qm2+tenth*qm2
      if(ik2(4) /= (id2(4)+iqmm2)) return
      call c0t5s(bd1(1),bd1(3),qm1,bk1(1),bk1(3),a2)
      if(dabs(a2) < eps) return
      call c0t5s(bd2(1),bd2(3),qm2,bk2(1),bk2(3),a3)
      if(dabs(a3) < eps) return
      call rmeajj(ik1(3),ik1(1),ik1(7),ik1(6),id1(1),id1(7),id1(6),s1)
      if(dabs(s1) < eps) return
      call rmeajj(ik2(3),ik2(1),ik2(7),ik2(6),id2(1),id2(7),id2(6),s2)
      if(dabs(s2) < eps) return
      a1=s1*s2*a2*a3
      call recop2(ns,ia,ib,ks1,ks2,ka,1,iat,rec)
      coeff=a1*rec/dsqrt(dble((2*ka+1)*(ik1(7)+1)*(ik2(7)+1)))
      jib=ib-1
      ifaz=0
      do  i=ia,jib
        ij=jlist(i)
        ifaz=ifaz+nq1(ij)
      end do
      ifaz=ifaz+1
      if(mod(ifaz,2) /= 0)coeff=-coeff
      if(ia.ne.ja) then
        ifaz = ik1(3)+ik2(3)-2*ka+2
        if(mod(ifaz,4) /= 0)coeff=-coeff
      endif
      coeff=-coeff*sqrt(dble(id1(3)+1))
      return
      end subroutine oneparticlejj2

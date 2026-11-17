module onescalar2_i
    implicit none (type, external)
    private
    public :: onescalar2
contains
!*******************************************************************
!                                                                  *
      subroutine onescalar2(jja,jjb,ja,jb,coeff)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
!                                                                  *
!     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
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
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, half, eps
      use m_c, only: double, jjc1, jjc2, jjq1, jjq2, jlist, klist, max_orbitals, ncore, npeel, nq1, nq2
      use trk_c, only: bd1, bd2, bd3, bd4, bk1, bk2, bk3, bk4, double, id1, id2, id3, id4, ik1, ik2, ik3, ik4
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use recoonescalar_i, only: recoonescalar
      use perko2_i, only: perko2
      use reco2_i, only: reco2
      use gg12_i, only: gg12
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb
      real(real64), intent(out) :: coeff
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: iat,jaa,jbb,nn,ib1,ii,inn,ifaz
      real(real64) :: qm1,qm2,rec,ww
!-----------------------------------------------
      coeff=zero
      if(ja == jb) return
      if(ja < jb) then
        jaa=ja
        jbb=jb
      else
        jaa=jb
        jbb=ja
      end if
      call recoonescalar(-1,jaa,jbb,jbb,jbb,1,iat)
      if(iat == 0)return
      qm1=half
      qm2=-half
      call perko2(ja,jb,ja,ja,2)
      if(id1(3) /= id2(3)) return
      call reco2(jaa,jbb,id2(3),0,iat,rec)
      if(iat == 0)return
      call gg12(ik1,ik2,bk1,bk2,id1,id2,bd1,bd2,qm1,qm2,ww)
      if(dabs(ww) > eps) then
         call reco2(jaa,jbb,id2(3),1,iat,rec)
         coeff=ww*rec*dsqrt(dble(id1(3)+1))
         nn=0
         ib1=jbb-1
         do ii=jaa,ib1
           inn=jlist(ii)
           nn=nq1(inn)+nn
         end do
         if((nn/2)*2 == nn) coeff=-coeff
!GG         IF(JA.GT.JB) COEFF=-COEFF
         coeff=-coeff
!
!     TRANSFORM FANO & RACAH PHASE CONVENTION
!     TO CONDON & SHORTLEY PHASE CONVENTION
!
        ifaz=ik1(5)*ik1(4)+ik2(5)*ik2(4)-id1(5)*id1(4)-id2(5)*id2(4)
        if((ifaz/4)*4 /= ifaz)coeff=-coeff
      end if
      return
      end subroutine onescalar2
end module onescalar2_i

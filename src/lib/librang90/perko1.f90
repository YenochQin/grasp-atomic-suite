!*******************************************************************
!                                                                  *
      subroutine perko1(ja,bk,ik,bd,id)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 22  -------------  *
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (FOR ONE SHELL)    *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
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
      use cons_c,          only: half
      use m_c,             only: jlist, nq1, nq2, jjq1, jjq2
      use orb_c,           only: orbital_principal_n, orbital_kappa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use nmtejj_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)                :: ja
      integer,      intent(out), dimension(7) :: ik, id
      real(double), intent(out), dimension(3) :: bk, bd
!      DIMENSION BK(3),IK(7),BD(3),ID(7)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij
!-----------------------------------------------
      ij=jlist(ja)
      ik(2)=orbital_principal_n(ij)
      id(2)=ik(2)
      ik(3)=(iabs(orbital_kappa(ij))*2)-1
      id(3)=ik(3)
      ik(4)=nq1(ij)
      id(4)=nq2(ij)
      ik(5)=(ik(3)+orbital_kappa(ij)/iabs(orbital_kappa(ij)))/2
      id(5)=ik(5)
      ik(6)=jjq1(3,ij)-1
      id(6)=jjq2(3,ij)-1
      ik(7)=iabs(orbital_kappa(ij))-jjq1(1,ij)
      id(7)=iabs(orbital_kappa(ij))-jjq2(1,ij)
      bk(1)=half*dble(ik(7))
      bd(1)=half*dble(id(7))
      bk(2)=half*dble(ik(6))
      bd(2)=half*dble(id(6))
      bk(3)=-half*dble(iabs(orbital_kappa(ij))-ik(4))
      bd(3)=-half*dble(iabs(orbital_kappa(ij))-id(4))
      ik(1)=nmtejj(ik(7),ik(6),ik(3),id(4),ik(4))
      id(1)=nmtejj(id(7),id(6),id(3),id(4),ik(4))
      return
      end subroutine perko1

!***********************************************************************
!                                                                      *
      subroutine consis(j)
!                                                                      *
!   This routine computes the weighted self-consistency of orbital J   *
!                                                                      *
!   Written by Farid A Parpia, at OXFORD    Last update: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use int_c
      use scf_c
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtp, i
      real(double) :: scmea, deltao
!-----------------------------------------------
!
!
      scmea = 0.0d00
      mtp = min(mtp0,mf(j))
      do i = 1, mtp
         deltao = abs(p(i)-pf(i,j)) + abs(q(i)-qf(i,j))
         scmea = dmax1(deltao,scmea)
      end do
      scnsty(j) = scmea*sqrt(ucf(j))
!
      return
      end subroutine consis

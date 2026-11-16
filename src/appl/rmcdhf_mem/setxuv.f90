!***********************************************************************
!                                                                      *
      subroutine setxuv(j)
!                                                                      *
!   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
!   subprograms IN and  OUT.                                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use grid_c
      use int_c
      use pote_c
      use scf_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nm1, i
      real(double) :: dmhh
!-----------------------------------------------
!
!   Define constants
!
      dmhh = -h*0.5d00
!
!   Set up arrays XU and XV; since XU(1), XV(1) are never used,
!   set them to some arbitrary value
!
      nm1 = n - 1
      xu(1) = 0.0d00
      xv(1) = 0.0d00
      xu(2:nm1) = dmhh*(xp(3:nm1+1)*rpor(3:nm1+1)+xp(2:nm1)*rpor(2:nm1)) + dp(2&
         :nm1)
      xv(2:nm1) = dmhh*(xq(3:nm1+1)*rpor(3:nm1+1)+xq(2:nm1)*rpor(2:nm1)) + dq(2&
         :nm1)
!
      xu(n) = dmhh*xp(n)*rpor(n) + dp(n)
      xv(n) = dmhh*xp(n)*rpor(n) + dq(n)
!
      return
      end subroutine setxuv

!***********************************************************************
!                                                                      *
      subroutine setxv(j)
!                                                                      *
!   This  subprogram  sets up the inhomogeneous terms for the varia-   *
!   tion equations.                                                    *
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
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: hhc
!-----------------------------------------------
!
      hhc = 0.5d00*h/c
!
!   Set up arrays TF and TG
!
      do i = 1, n
         xu(i) = -qf(i,j)*hhc*rp(i)
         xv(i) = pf(i,j)*hhc*rp(i)
      end do
!
      xu(:n-1) = xu(2:n) + xu(:n-1)
      xv(:n-1) = xv(2:n) + xv(:n-1)
!
      return
!
      end subroutine setxv

!***********************************************************************
!                                                                      *
        subroutine genintbreit1wrap (myid, nprocs, j2max)

!   Written by Per Jonsson                Last revision: October 2014  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use bilst_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintbreit1_i
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: gisummpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: myid, nprocs
      integer             :: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n
!-----------------------------------------------
!
      call genintbreit1 ((myid), (nprocs), n, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      call gisummpi (indtp1, n)
      call gdsummpi (valtp1, n)

      return
      end subroutine genintbreit1wrap

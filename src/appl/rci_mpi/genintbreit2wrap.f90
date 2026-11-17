!***********************************************************************
!                                                                      *
        subroutine genintbreit2wrap (myid, nprocs, j2max)

!   Written by Per Jonsson                Last revision: October 2014  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use bilst_c, only: double, ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintbreit2_i, only: genintbreit2
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: gisummpi
      implicit none (type, external)
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
      call genintbreit2 ((myid), (nprocs), n, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      call gisummpi (indtp2, n)
      call gdsummpi (valtp2, n)

      return
      end subroutine genintbreit2wrap

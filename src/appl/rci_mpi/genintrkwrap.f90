!***********************************************************************
!                                                                      *
      subroutine genintrkwrap(myid, nprocs, j2max)
!
!   Written by     Xinghong He            Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:42:46   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cteilsrk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintrk_i
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: gisummpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: myid
      integer  :: nprocs
      integer  :: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: n
!-----------------------------------------------

      call genintrk (myid, nprocs, n, j2max)

! Gather integrals (and their indeces) from- and send to- all nodes

      call gisummpi (indteirk, n)
      call gdsummpi (valteirk, n)

      return
      end subroutine genintrkwrap

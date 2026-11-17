!***********************************************************************
!                                                                      *
      subroutine genintrkwrap(myid, nprocs, j2max)
!
!   Written by     Xinghong He            Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use cteilsrk_c, only: double, indteirk, valteirk, rk_local_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintrk_i, only: genintrk
      use mpiu_i, only: gdsummpi
      use mpiu_i, only: gisummpi
      implicit none (type, external)
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

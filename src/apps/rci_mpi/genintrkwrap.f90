!***********************************************************************
!                                                                      *
      subroutine genintrkwrap(mpi_rank, mpi_process_count, max_orbital_two_j)
!
!   Written by     Xinghong He            Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintrk_i, only: genintrk
      use mpi_utilities, only: gdsummpi
      use mpi_utilities, only: gisummpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: mpi_rank
      integer  :: mpi_process_count
      integer  :: max_orbital_two_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: n
!-----------------------------------------------

      call genintrk (mpi_rank, mpi_process_count, n, max_orbital_two_j)

! Gather integrals (and their indeces) from- and send to- all nodes

      call gisummpi (indteirk, n)
      call gdsummpi (valteirk, n)

      return
      end subroutine genintrkwrap

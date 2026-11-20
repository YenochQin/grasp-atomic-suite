!***********************************************************************
!                                                                      *
        subroutine genintbreit2wrap (mpi_rank, mpi_process_count, max_orbital_two_j)

!   Written by Per Jonsson                Last revision: October 2014  *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use genintbreit2_i, only: genintbreit2
      use mpi_utilities, only: gdsummpi
      use mpi_utilities, only: gisummpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mpi_rank, mpi_process_count
      integer             :: max_orbital_two_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n
!-----------------------------------------------
!
      call genintbreit2 ((mpi_rank), (mpi_process_count), n, max_orbital_two_j)

! Gather integrals (and their indeces) from- and send to- all nodes

      call gisummpi (indtp2, n)
      call gdsummpi (valtp2, n)

      return
      end subroutine genintbreit2wrap

      subroutine cofpot(extended_optimal_level, orbital_index, radial_point_count)
!-----------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use scf_radial_potentials, only: max_radial_points, xp, xq, yp
      use mpi_s, only: mpi_process_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i, only: setcof
      use ypot_i, only: ypot
      use xpot_i, only: xpot
      use lagcon_i, only: lagcon
      use dacon_i, only: dacon
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: orbital_index
      integer  :: radial_point_count
      logical  :: extended_optimal_level
!-----------------------------------------------------------------------

      call setcof (extended_optimal_level, orbital_index)
      call ypot (orbital_index)
      call xpot (orbital_index)
      call lagcon (orbital_index, mpi_process_count)
      call dacon

      return
      end subroutine cofpot

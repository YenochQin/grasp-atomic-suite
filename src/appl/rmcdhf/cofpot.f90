      subroutine cofpot(eol, j, npts)
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
      integer  :: j
      integer  :: npts
      logical  :: eol
!-----------------------------------------------------------------------

      call setcof (eol, j)
      call ypot (j)
      call xpot (j)
      call lagcon (j, mpi_process_count)
      call dacon

      return
      end subroutine cofpot

      subroutine cofpotmpi(eol, j, npts)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use pote_c, only: double, max_radial_points, xp, xq, yp
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i, only: setcof
      use ypot_i, only: ypot
      use xpot_i, only: xpot
      use lagcon_i, only: lagcon
      use dacon_i, only: dacon
      use pote_c, only: double, max_radial_points, xp, xq, yp
      use suite_linear_algebra, only: dcopy
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
      integer  :: npts
      logical  :: eol
      real(real64) :: tmpmpi(npts)
!-----------------------------------------------------------------------

      call setcof (eol, j)
      call ypot (j)
      call xpot (j)
      call lagcon (j, mpi_process_count)
      call dacon

      call mpi_allreduce (yp, tmpmpi, npts, mpi_double_precision,  &
                          mpi_sum, mpi_comm_world, mpi_error_code)
      call dcopy (npts, tmpmpi, 1, yp, 1)

      call mpi_allreduce (xp, tmpmpi, npts, mpi_double_precision,  &
                          mpi_sum, mpi_comm_world, mpi_error_code)
      call dcopy (npts, tmpmpi, 1, xp, 1)

      call mpi_allreduce (xq, tmpmpi, npts, mpi_double_precision,  &
                          mpi_sum, mpi_comm_world, mpi_error_code)
      call dcopy (npts, tmpmpi, 1, xq, 1)

      return
      end subroutine cofpotmpi

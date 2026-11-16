      subroutine cofpotmpi(eol, j, npts)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use pote_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i
      use ypot_i
      use xpot_i
      use lagcon_i
      use dacon_i
      use pote_c
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

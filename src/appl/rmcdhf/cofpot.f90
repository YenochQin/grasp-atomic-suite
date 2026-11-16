      subroutine cofpot(eol, j, npts)
!-----------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use pote_c
      use mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i
      use ypot_i
      use xpot_i
      use lagcon_i
      use dacon_i
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

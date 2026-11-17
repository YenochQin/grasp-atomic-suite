!***********************************************************************
      subroutine setmcp2(myid, nprocs, ncore, idblk, filehead)
!
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
! and MPI programs.
!
!   Written by Xinghong He                Last revision: 30 Jun 1998
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use def_c, only: double, electron_count
      use iccu_c, only: iccut, iccutblk
      use mcp_c, only: diagonal_only, double, max_coulomb_multipole_rank, use_first_order_restriction
      use orb_c, only: csf_count, double, orbital_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setmcp_i, only: setmcp
      use getinf_i, only: getinf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: myid
      integer  :: nprocs
      integer  :: ncore
      character  :: filehead*(*)
      character  :: idblk(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k

      call setmcp (myid, nprocs, ncore, idblk, filehead)

! DIAG, ICCUTblk, LFORDR are set in GETINF

      if (myid == 0) call getinf

      do k = 30, 32 + max_coulomb_multipole_rank
         write (k) electron_count, csf_count, orbital_count
         write (k) diagonal_only, iccut(1), use_first_order_restriction
      end do

      do k = 30, 32 + max_coulomb_multipole_rank
         close(k)
      end do

      return
      end subroutine setmcp2

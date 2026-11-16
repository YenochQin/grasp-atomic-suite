!***********************************************************************
      subroutine setmcp2(myid, nprocs, ncore, idblk, filehead)
!
! A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
! and MPI programs.
!
!   Written by Xinghong He                Last revision: 30 Jun 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      use def_c
      use iccu_c
      use mcp_c
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setmcp_i
      use getinf_i
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: myid
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

      do k = 30, 32 + kmax
         write (k) nelec, ncf, nw
         write (k) diag, iccut(1), lfordr
      end do

      do k = 30, 32 + kmax
         close(k)
      end do

      return
      end subroutine setmcp2

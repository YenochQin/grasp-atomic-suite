!***********************************************************************
!                                                                      *
      subroutine setdbgmpi(dbgfile)
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 10 Dec 1992   *
!   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setdbg_i, only: setdbg
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: dbgfile
!-----------------------------------------------
!
      if (mpi_rank == 0) then
         call setdbg (dbgfile)
      end if

      call mpi_bcast (ldbpa, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast (ldbpg, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast (ldbpr,30, mpi_logical, 0, mpi_comm_world, mpi_error_code)

      return
      end subroutine setdbgmpi

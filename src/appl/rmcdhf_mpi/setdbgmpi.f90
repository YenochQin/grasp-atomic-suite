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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:22:29   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setdbg_i, only: setdbg
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: dbgfile
!-----------------------------------------------
!
      if (mpi_rank .eq. 0) then
         call setdbg (dbgfile)
      endif

      call mpi_bcast (ldbpa, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast (ldbpg, 5, mpi_logical, 0, mpi_comm_world, mpi_error_code)
      call mpi_bcast (ldbpr,30, mpi_logical, 0, mpi_comm_world, mpi_error_code)

      return
      end subroutine setdbgmpi

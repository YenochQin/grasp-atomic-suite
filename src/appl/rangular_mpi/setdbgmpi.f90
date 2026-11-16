!***********************************************************************
!                                                                      *
      subroutine setdbgmpi (debug, fullname)
!                                                                      *
!   This routine calls setdbg to open the  .dbg  file and sets the     *
!   arrays that control debug printout from the GENMCP program.        *
!                                                                      *
!   Written by Xinghong He                  Last update: 29 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setdbg_i
      use debug_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: debug
      character(len=*) :: fullname
!-----------------------------------------------
!
     if (mpi_rank == 0) then
         call setdbg (debug, fullname)
      endif
      call mpi_bcast (debug,1,mpi_logical,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ldbpa,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ldbpg,5,mpi_logical,0,mpi_comm_world,mpi_error_code)
      return
      end subroutine setdbgmpi

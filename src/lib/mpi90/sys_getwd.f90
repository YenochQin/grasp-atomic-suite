!***********************************************************************
      subroutine sys_getwd (dir,lcwd)
!
!      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!
!   Xinghong He  98-08-24
!
!************************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_C
      use suite_system, only: working_directory

      USE mpiu_I, ONLY: stopmpi
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: lcwd

      INTEGER status


      status = working_directory(dir)
      lcwd = len_trim(dir)

      if (status.ne.0) then
!        print*, 'couldn''t get the current directory, exiting...'
!cjb myid =
         CALL stopmpi ('sys_getwd: status.ne.0; myid = ', myid)

      end if
      return
      end subroutine sys_getwd

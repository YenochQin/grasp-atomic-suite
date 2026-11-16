module sys_getwd_i
    implicit none (type, external)
    private
    public :: sys_getwd
contains
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c
      use suite_system, only: working_directory

      use mpi_runtime_i, only: stopmpi
      implicit none (type, external)
      character(len=*), intent(out):: dir
      integer, intent(out):: lcwd

      integer :: status


      status = working_directory(dir)
      lcwd = len_trim(dir)

      if (status/=0) then
!        print*, 'couldn''t get the current directory, exiting...'
!cjb myid =
         call stopmpi ("sys_getwd: status.ne.0; myid = ", mpi_rank)

      end if
      return
      end subroutine sys_getwd
end module sys_getwd_i

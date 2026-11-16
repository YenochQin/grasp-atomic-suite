module mpiu_i
    implicit none (type, external)
    interface
      subroutine startmpi (myid, nprocs, host, lenhost)
      implicit none (type, external)
      integer, intent(out)  :: myid,  nprocs, lenhost
      character(len=*), intent(out) :: host
      end subroutine startmpi
      subroutine mpix_startup (myid, nprocs, host, lenhost, ncount1, &
                               progname, description, infiledesc, &
                               outfiledesc, quiet_workers)
      implicit none (type, external)
      integer :: myid, nprocs, lenhost, ncount1
      character(len=*), intent(out) :: host
      character(len=*), intent(in) :: progname, description
      character(len=*), intent(in) :: infiledesc, outfiledesc
      logical, intent(in) :: quiet_workers
      end subroutine mpix_startup
      subroutine startmpi2 (myid, nprocs, host, lenhost, ncount1, &
                           startdir, permdir, tmpdir, progname)
      implicit none (type, external)
      integer, intent(out)  :: myid,  nprocs, lenhost
      integer, intent(out)  :: ncount1
      character(len=*), intent(out) :: host
      character(len=*), intent(in) :: progname
      character(len=*), intent(out) :: startdir, permdir, tmpdir
      end subroutine startmpi2
      subroutine mpix_shutdown (myid, ncount1, progname)
      implicit none (type, external)
      integer, intent(in) :: myid, ncount1
      character(len=*), intent(in) :: progname
      end subroutine mpix_shutdown
      subroutine stopmpi (what, myid)
      implicit none (type, external)
      integer, intent(in)  ::  myid
      character(len=*), intent(in) :: what
      end subroutine stopmpi
      subroutine stopmpi2 (myid, nprocs, host, lenhost, ncount1, &
                            progname)
      implicit none (type, external)
      integer, intent(in)  ::  myid,nprocs,lenhost
      integer              ::  ncount1
      character*(*), intent(in) ::  host, progname
      end subroutine stopmpi2
      subroutine mpix_printmsg (msg, myid, nprocs)
      implicit none (type, external)
      integer, intent(in)  ::  myid, nprocs
      character(len=*), intent(inout) :: msg
      end subroutine mpix_printmsg
      subroutine mpix_chkpt (myid, what)
      implicit none (type, external)
      integer, intent(in) :: myid
      character(len=*), intent(in) :: what
      end subroutine mpix_chkpt
      subroutine mpix_bytes (n, newtype, ierr)
      implicit none (type, external)
      integer, intent(in) :: n
      integer, intent(out) :: newtype, ierr
      end subroutine mpix_bytes
      subroutine gdsummpi (x, n)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)                         :: n
      real(real64), intent(inout) ::   x(1:n)
      end subroutine gdsummpi
      subroutine gdrsummpi_root (x, n)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)                          :: n
      real(real64), intent(inout) :: x(1:n)
      end subroutine gdrsummpi_root
      subroutine gdmaxmpi_root (x, n)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in)                          :: n
      real(real64), intent(inout) :: x(1:n)
      end subroutine gdmaxmpi_root
      subroutine gisummpi (ix, n)
      implicit none (type, external)
      integer, intent(in)                         :: n
      integer, intent(inout) ::  ix(1:n)
      end subroutine gisummpi
    end interface
end module

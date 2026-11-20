      module mcpmpi_i
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine mcpmpi (nb, is_restart, output_file_prefix)
      implicit none (type, external)
      integer  :: nb
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      logical, intent(in) :: is_restart
      character(len=*), intent(in) :: output_file_prefix
      end subroutine mcpmpi
      end interface
      end module mcpmpi_i

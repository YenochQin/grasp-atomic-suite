      module mcpmpi_i
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine mcpmpi (nb, restrt, fhead)
      implicit none (type, external)
      integer  :: nb
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      logical, intent(in) :: restrt
      character(len=*), intent(in) :: fhead
      end subroutine mcpmpi
      end interface
      end module mcpmpi_i

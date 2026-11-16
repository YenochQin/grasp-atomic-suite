      module mcpmpi_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!cjb  myid, nprocs = NOT args
      subroutine mcpmpi (nb, restrt, fhead)
      implicit none (type, external)
      integer  :: nb
!     INTEGER  :: MYID
!     INTEGER  :: NPROCS
      logical, intent(in) :: restrt
      character(len=*), intent(in) :: fhead
      end subroutine
      end interface
      end module

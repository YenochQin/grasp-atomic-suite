      module setmcpmpi_i
      implicit none (type, external)
      interface
!cjb  myid, nprocs = NOT args
      subroutine setmcpmpi (ncore, idblk, filehead)
      implicit none (type, external)
!     INTEGER , INTENT(IN) :: MYID
!     INTEGER , INTENT(IN) :: NPROCS
      integer , intent(in) :: ncore
      character (len = 8), intent(in) :: idblk(*)
      character (len = *), intent(in) :: filehead
!...This routine performs I/O.
      end subroutine setmcpmpi
      end interface
      end module setmcpmpi_i

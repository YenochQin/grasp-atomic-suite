      module mcp_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine mcp (nb, restrt, myid, nprocs, fhead)
      integer  :: nb
      integer  :: myid
      integer  :: nprocs
      logical, intent(in) :: restrt
      character(len=*), intent(in) :: fhead
      end subroutine
      end interface
      end module

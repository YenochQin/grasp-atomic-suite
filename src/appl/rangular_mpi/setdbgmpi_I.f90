      module setdbgmpi_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine setdbgmpi (debug, fullname)
      implicit none (type, external)
      logical :: debug
      character(len=*) :: fullname
      end subroutine setdbgmpi
      end interface
      end module

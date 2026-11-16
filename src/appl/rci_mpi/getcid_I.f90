      module getcid_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine getcid (isofile, rwffile, idblk)
      implicit none (type, external)
      character(len=*):: isofile, rwffile
      character(len=8) ::  idblk(*)
      end subroutine
      end interface
      end module

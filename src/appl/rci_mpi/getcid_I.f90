      module getcid_i
      implicit none (type, external)
      interface
      subroutine getcid (isofile, rwffile, idblk)
      implicit none (type, external)
      character(len=*):: isofile, rwffile
      character(len=8) ::  idblk(*)
      end subroutine getcid
      end interface
      end module getcid_i

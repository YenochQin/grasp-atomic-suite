      module setres_i
      implicit none (type, external)
      interface
      subroutine setres (isofile, rwffile, idblk)
      implicit none (type, external)
      character (len = *) :: isofile
      character (len = *) :: rwffile
      character (len = 8) :: idblk(*)
      end subroutine
      end interface
      end module

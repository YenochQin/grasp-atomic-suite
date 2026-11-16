      module getscdmpi_i
      implicit none (type, external)
      interface
      subroutine getscdmpi (eol, idblk, isofile, rwffile)
      implicit none (type, external)
      logical, intent(out) :: eol
      character (len = 8) :: idblk(*)
      character (len = *) :: isofile
      character (len = *) :: rwffile
      end subroutine getscdmpi
      end interface
      end module getscdmpi_i

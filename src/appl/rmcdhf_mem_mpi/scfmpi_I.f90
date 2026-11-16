      module scfmpi_i
      implicit none (type, external)
      interface
      subroutine scfmpi (eol, rwffile2)
      implicit none (type, external)
      logical, intent(in) :: eol
      character (len = *) :: rwffile2
      end subroutine
      end interface
      end module

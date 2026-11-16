      module scf_i
      implicit none (type, external)
      interface
      subroutine scf (eol, rwffile2)
      implicit none (type, external)
      logical, intent(in) :: eol
      character (len = *) :: rwffile2
      end subroutine scf
      end interface
      end module scf_i

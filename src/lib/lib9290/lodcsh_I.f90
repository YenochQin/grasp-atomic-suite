      module lodcsh_i
      implicit none (type, external)
      interface
      subroutine lodcsh (nfile, ncore)
      implicit none (type, external)
      integer, intent(in) :: nfile
      integer, intent(out) :: ncore
      end subroutine
      end interface
      end module

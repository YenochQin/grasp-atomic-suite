      module lodcsh2gg_i
      implicit none (type, external)
      interface
      subroutine lodcsh2gg(nfile, ncore, jb)
      implicit none (type, external)
      integer, intent(in) :: nfile
      integer, intent(out) :: ncore
      integer , intent(in) :: jb
      end subroutine
      end interface
      end module

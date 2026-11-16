      module lodcsh2_i
      implicit none (type, external)
      interface
      subroutine lodcsh2(nfile, ncore, jb)
      implicit none (type, external)
      integer, intent(in) :: nfile
      integer, intent(in) :: ncore
      integer , intent(in) :: jb
      end subroutine
      end interface
      end module

      module mctout_i
      implicit none (type, external)
      interface
      subroutine mctout (iopar, jkp, name)
      implicit none (type, external)
      integer, intent(out) :: iopar
      integer, intent(in) :: jkp
      character (len = 24) :: name(2)
      end subroutine
      end interface
      end module

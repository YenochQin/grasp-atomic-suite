      module merg12_i
      implicit none (type, external)
      interface
      subroutine merg12 (name, ncorer, ncore)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
      integer, intent(in) :: ncorer
      integer, intent(in) :: ncore
      end subroutine
      end interface
      end module

      module merg12_i
      implicit none (type, external)
      interface
      subroutine merg12 (name, ncorer, ncore)
      implicit none (type, external)
      character (len = 128), intent(in) :: name(2)
      integer, intent(in) :: ncorer
      integer, intent(in) :: ncore
      end subroutine merg12
      end interface
      end module merg12_i

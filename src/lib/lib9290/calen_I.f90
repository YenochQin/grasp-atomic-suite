      module calen_i
      implicit none (type, external)
      interface
      subroutine calen (jtime, jdate)
      implicit none (type, external)
      character(len=10), intent(out) :: jtime
      character(len=8), intent(out) :: jdate
      end subroutine
      end interface
      end module

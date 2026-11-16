      module outsdampi_i
      implicit none (type, external)
      interface
      subroutine outsdampi (lprint, nnonz, ncf)
      implicit none (type, external)
      logical, intent(in) :: lprint
      integer, intent(in) :: nnonz
      integer, intent(in) :: ncf
      end subroutine
      end interface
      end module

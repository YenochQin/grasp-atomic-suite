      module tcsl_i
      implicit none (type, external)
      interface
      subroutine tcsl (n)
      implicit none (type, external)
      integer, intent(in) :: n
      end subroutine
      end interface
      end module

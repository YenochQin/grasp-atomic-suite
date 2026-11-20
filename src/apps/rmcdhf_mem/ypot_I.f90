      module ypot_i
      implicit none (type, external)
      interface
      subroutine ypot (orbital_index)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      end subroutine ypot
      end interface
      end module ypot_i

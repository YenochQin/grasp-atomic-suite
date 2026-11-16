      module ypot_i
      implicit none (type, external)
      interface
      subroutine ypot (j)
      implicit none (type, external)
      integer, intent(in) :: j
      end subroutine ypot
      end interface
      end module ypot_i

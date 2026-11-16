      module packls_i
      implicit none (type, external)
      interface
      subroutine packls (m, el, q, couple, str)
      implicit none (type, external)
      integer, intent(in) :: m
      character (len = 3), intent(in) :: el(*)
      integer, intent(in) :: q(*)
      character (len = 3), intent(in) :: couple(*)
      character (len = 164), intent(out) :: str
      end subroutine packls
      end interface
      end module packls_i

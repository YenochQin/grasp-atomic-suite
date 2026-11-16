      module jthn_i
      implicit none (type, external)
      interface
!
      integer function jthn(k,n,i)
      implicit none (type, external)
      integer, intent(in) :: i, n, k
      end function
      end interface
      end module

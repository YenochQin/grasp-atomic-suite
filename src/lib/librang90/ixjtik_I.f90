      module ixjtik_i
      implicit none (type, external)
      interface
      integer function ixjtik (i, j, k, l, m, n)
      implicit none (type, external)
      integer :: i
      integer :: j
      integer :: k
      integer :: l
      integer :: m
      integer :: n
!VAST...Calls: ITTK
      end function
      end interface
      end module

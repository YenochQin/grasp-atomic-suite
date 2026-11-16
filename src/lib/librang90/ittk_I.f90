      module ittk_i
      implicit none (type, external)
      interface
      integer function ittk (i, j, k)
      implicit none (type, external)
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer, intent(in) :: k
      end function
      end interface
      end module

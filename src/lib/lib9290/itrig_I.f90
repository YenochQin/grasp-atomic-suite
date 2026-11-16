      module itrig_i
      implicit none (type, external)
      interface
      integer function itrig (i1, i2, i3)
      implicit none (type, external)
      integer, intent(in) :: i1
      integer, intent(in) :: i2
      integer, intent(in) :: i3
      end function
      end interface
      end module

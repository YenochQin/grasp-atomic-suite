      module jfaze_i
      implicit none (type, external)
      interface
!
      integer function jfaze(i1,i2,i3,i4)
      implicit none (type, external)
      integer, intent(in) :: i1, i2, i3, i4
      end function
      end interface
      end module

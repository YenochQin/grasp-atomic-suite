      module itrexg_i
      implicit none (type, external)
      interface
!
      integer function itrexg(i1,i2,i3,i4,k)
      implicit none (type, external)
      integer, intent(in)  :: i1, i2, i3, i4
      integer, intent(out) :: k
      end function
      end interface
      end module

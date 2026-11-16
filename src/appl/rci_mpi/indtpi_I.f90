      module indtpi_i
      implicit none (type, external)
      interface
      integer function indtpi (itype, i)
      implicit none (type, external)
      integer, intent(in) :: itype
      integer :: i
      end function
      end interface
      end module

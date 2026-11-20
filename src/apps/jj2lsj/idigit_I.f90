      module idigit_i
      implicit none (type, external)
      interface
      integer function idigit (cst)
      implicit none (type, external)
      character (len = 1), intent(in) :: cst
      end function idigit
      end interface
      end module idigit_i

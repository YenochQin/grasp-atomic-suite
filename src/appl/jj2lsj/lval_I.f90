      module lval_i
      implicit none (type, external)
      interface
      integer function lval (symbol)
      implicit none (type, external)
      character (len = 1), intent(in) :: symbol
      end function
      end interface
      end module

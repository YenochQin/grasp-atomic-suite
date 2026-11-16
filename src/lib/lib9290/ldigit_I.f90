      module ldigit_i
      implicit none (type, external)
      interface
      logical function ldigit (cst)
      implicit none (type, external)
      character (len = 1), intent(in) :: cst
      end function
      end interface
      end module

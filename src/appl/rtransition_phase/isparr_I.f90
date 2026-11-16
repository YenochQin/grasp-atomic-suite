      module isparr_i
      implicit none (type, external)
      interface
      integer function isparr (icsf)
      implicit none (type, external)
      integer, intent(in) :: icsf
      end function
      end interface
      end module

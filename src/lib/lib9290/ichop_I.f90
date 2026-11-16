      module ichop_i
      implicit none (type, external)
      interface
      integer function ichop (isubsh, icsf)
      implicit none (type, external)
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
      end function
      end interface
      end module

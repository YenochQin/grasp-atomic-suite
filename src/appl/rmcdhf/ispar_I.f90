      module ispar_i
      implicit none (type, external)
      interface
      integer function ispar (icsf)
      implicit none (type, external)
      integer, intent(in) :: icsf
      end function ispar
      end interface
      end module ispar_i

      module iqr_i
      implicit none (type, external)
      interface
      integer function iqr (isubsh, icsf)
      implicit none (type, external)
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
      end function iqr
      end interface
      end module iqr_i

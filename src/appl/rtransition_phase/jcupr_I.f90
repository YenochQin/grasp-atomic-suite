      module jcupr_i
      implicit none (type, external)
      interface
      integer function jcupr (loc, icsf)
      implicit none (type, external)
      integer, intent(in) :: loc
      integer, intent(in) :: icsf
      end function jcupr
      end interface
      end module jcupr_i

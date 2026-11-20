      module jcupr_i
      implicit none (type, external)
      interface
      integer function jcupr (loc, csf_index)
      implicit none (type, external)
      integer, intent(in) :: loc
      integer, intent(in) :: csf_index
      end function jcupr
      end interface
      end module jcupr_i

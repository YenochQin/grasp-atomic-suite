      module isparr_i
      implicit none (type, external)
      interface
      integer function isparr (csf_index)
      implicit none (type, external)
      integer, intent(in) :: csf_index
      end function isparr
      end interface
      end module isparr_i

      module iqr_i
      implicit none (type, external)
      interface
      integer function iqr (subshell_index, csf_index)
      implicit none (type, external)
      integer, intent(in) :: subshell_index
      integer, intent(in) :: csf_index
      end function iqr
      end interface
      end module iqr_i

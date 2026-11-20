      module jqsr_i
      implicit none (type, external)
      interface
      integer function jqsr (iwhich, subshell_index, csf_index)
      implicit none (type, external)
      integer, intent(in) :: iwhich
      integer, intent(in) :: subshell_index
      integer, intent(in) :: csf_index
      end function jqsr
      end interface
      end module jqsr_i

      module jqsr_i
      implicit none (type, external)
      interface
      integer function jqsr (iwhich, isubsh, icsf)
      implicit none (type, external)
      integer, intent(in) :: iwhich
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
      end function
      end interface
      end module

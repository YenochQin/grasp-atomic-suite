      module jqs_i
      implicit none (type, external)
      interface
      integer function jqs (iwhich, isubsh, icsf)
      implicit none (type, external)
      integer, intent(in) :: iwhich
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
!VAST...Calls: JQSA
      end function
      end interface
      end module

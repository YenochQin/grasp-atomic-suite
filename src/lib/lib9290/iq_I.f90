      module iq_i
      implicit none (type, external)
      interface
      integer function iq (isubsh, icsf)
      implicit none (type, external)
      integer, intent(in) :: isubsh
      integer, intent(in) :: icsf
!VAST...Calls: IQA
      end function
      end interface
      end module

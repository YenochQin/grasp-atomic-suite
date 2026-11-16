      module itjpor_i
      implicit none (type, external)
      interface
      integer function itjpor (icsf)
      implicit none (type, external)
      integer, intent(in) :: icsf
      end function
      end interface
      end module

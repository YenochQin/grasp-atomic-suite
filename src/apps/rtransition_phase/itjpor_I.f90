      module itjpor_i
      implicit none (type, external)
      interface
      integer function itjpor (csf_index)
      implicit none (type, external)
      integer, intent(in) :: csf_index
      end function itjpor
      end interface
      end module itjpor_i

      module itjpo_i
      implicit none (type, external)
      interface
      integer function itjpo (csf_index)
      implicit none (type, external)
      integer, intent(in) :: csf_index
      end function itjpo
      end interface
      end module itjpo_i

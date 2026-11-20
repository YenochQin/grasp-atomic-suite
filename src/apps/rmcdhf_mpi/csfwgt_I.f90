      module csfwgt_i
      implicit none (type, external)
      interface
      subroutine csfwgt (write_standard_output)
      implicit none (type, external)
      logical, intent(in) :: write_standard_output
      end subroutine csfwgt
      end interface
      end module csfwgt_i

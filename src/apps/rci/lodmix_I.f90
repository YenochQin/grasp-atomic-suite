      module lodmix_i
      implicit none (type, external)
      interface
      subroutine lodmix (block_labels)
      implicit none (type, external)
      character(len=8), intent(in):: block_labels(*)
      end subroutine lodmix
      end interface
      end module lodmix_i

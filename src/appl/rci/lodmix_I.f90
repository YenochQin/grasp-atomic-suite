      module lodmix_i
      implicit none (type, external)
      interface
      subroutine lodmix (idblk)
      implicit none (type, external)
      character(len=8), intent(in):: idblk(*)
      end subroutine lodmix
      end interface
      end module lodmix_i

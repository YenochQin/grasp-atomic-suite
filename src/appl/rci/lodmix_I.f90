      module lodmix_i
      implicit none (type, external)
      interface
      subroutine lodmix (idblk)
      implicit none (type, external)
      character(len=8), intent(in):: idblk(*)
      end subroutine
      end interface
      end module

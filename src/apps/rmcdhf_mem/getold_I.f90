      module getold_i
      implicit none (type, external)
      interface
      subroutine getold (idblk)
      implicit none (type, external)
      character (len = 8) :: idblk(*)
      end subroutine getold
      end interface
      end module getold_i

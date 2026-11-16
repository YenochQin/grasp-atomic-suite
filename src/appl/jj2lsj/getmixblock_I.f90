      module getmixblock_i
      implicit none (type, external)
      interface
      subroutine getmixblock (name, nci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(in) :: nci
      end subroutine getmixblock
      end interface
      end module getmixblock_i

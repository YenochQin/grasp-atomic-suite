      module getmix_i
      implicit none (type, external)
      interface
      subroutine getmix (name, inpci, iblk)
      implicit none (type, external)
      character (len = 128), intent(in) :: name
      integer, intent(in) :: inpci
      integer, intent(in) :: iblk
      end subroutine
      end interface
      end module

      module readmix_i
      implicit none (type, external)
      interface
      subroutine readmix (name, inpci, init)
      implicit none (type, external)
      character (len = 24) :: name(2)
      integer :: inpci
      integer, intent(in) :: init
      end subroutine
      end interface
      end module

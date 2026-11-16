      module mcpout_i
      implicit none (type, external)
      interface
      subroutine mcpout (name, startdir,ik, ntestg, inpci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      character (len = 128) :: startdir
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer :: inpci
      end subroutine
      end interface
      end module

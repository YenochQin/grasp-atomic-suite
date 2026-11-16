      module mcpout_i
      implicit none (type, external)
      interface
      subroutine mcpout (name, ik, ntestg, inpci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer, intent(in) :: ik
      integer, intent(in) :: ntestg
      integer :: inpci
      end subroutine mcpout
      end interface
      end module mcpout_i

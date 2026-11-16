      module genmcp_i
      implicit none (type, external)
      interface
      subroutine genmcp (name, ic, ntestg, inpci)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer :: ic
      integer :: ntestg
      integer :: inpci
      end subroutine
      end interface
      end module

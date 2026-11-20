      module setcslb_i
      implicit none (type, external)
      interface
      subroutine setcslb (name, core_subshell_count,igg)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer :: core_subshell_count,igg
      end subroutine setcslb
      end interface
      end module setcslb_i

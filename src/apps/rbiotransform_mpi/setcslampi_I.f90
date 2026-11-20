      module setcsla_i
      implicit none (type, external)
      interface
      subroutine setcsla (name, core_subshell_count,igg)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      integer :: core_subshell_count,igg
      end subroutine setcsla
      end interface
      end module setcsla_i

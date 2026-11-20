      module merg12_i
      implicit none (type, external)
      interface
      subroutine merg12 (state_name, initial_core_subshell_count, final_core_subshell_count)
      implicit none (type, external)
      character (len = 128), intent(in) :: state_name(2)
      integer, intent(in) :: initial_core_subshell_count
      integer, intent(in) :: final_core_subshell_count
      end subroutine merg12
      end interface
      end module merg12_i

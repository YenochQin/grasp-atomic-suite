      module getmixblock_i
      implicit none (type, external)
      interface
      subroutine getmixblock (state_name, ci_input_mode)
      implicit none (type, external)
      character (len = 24), intent(in) :: state_name
      integer, intent(in) :: ci_input_mode
      end subroutine getmixblock
      end interface
      end module getmixblock_i

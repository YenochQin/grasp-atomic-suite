      module ldcsl2_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine ldcsl2 (final_core_subshell_count, state_name)
      implicit none (type, external)
      integer, intent(out) :: final_core_subshell_count
      character (len = 128), intent(in) :: state_name
      end subroutine ldcsl2
      end interface
      end module ldcsl2_i

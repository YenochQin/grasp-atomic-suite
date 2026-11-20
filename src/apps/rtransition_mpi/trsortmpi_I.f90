      module trsort_i
      implicit none (type, external)
      interface
      subroutine trsort (name, nfile, nfile2, print_details, transition_type_index, iblki, iblkf)
      implicit none (type, external)
      character (len = 24), intent(in) :: name(2)
      integer, intent(in) :: nfile
      integer, intent(in) :: nfile2
      logical, intent(in) :: print_details
      integer :: transition_type_index
      integer, intent(in) :: iblki
      integer, intent(in) :: iblkf
      end subroutine trsort
      end interface
      end module trsort_i

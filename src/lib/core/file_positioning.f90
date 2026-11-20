! Position a Fortran unit by rewinding and advancing or by backspacing records.
module file_positioning
    implicit none (type, external)
    private
    public :: posfile
contains
      subroutine posfile(position_mode, file_unit, record_count)
      implicit none (type, external)
      integer, intent(in) :: position_mode
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_count
      integer, parameter :: rewind_and_advance = 0
      integer, parameter :: backspace_records = 1
      integer :: record_index

      select case (position_mode)

      case (rewind_and_advance)
         rewind (file_unit)
         do record_index = 1, record_count
            read (file_unit)
         end do

      case (backspace_records)
         do record_index = 1, record_count
            backspace (file_unit)
         end do

      end select

      return
      end subroutine posfile
end module file_positioning

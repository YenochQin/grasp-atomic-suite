! Position a Fortran unit by rewinding and advancing or by backspacing records.
module file_positioning
    implicit none (type, external)

    ! Module interface
    private
    public :: position_calculation_file
contains
    subroutine position_calculation_file(position_mode, file_unit, record_count)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: position_mode
        integer, intent(in) :: file_unit
        integer, intent(in) :: record_count

        ! Local declarations
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
    end subroutine position_calculation_file
end module file_positioning

! Check the isotope-file header on unit 22, load its records, and close the file.
module isotope_data_input
    implicit none (type, external)
    private
    public :: setiso
contains
      subroutine setiso(isotope_file_path)
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use file_open, only: openfl
      use isotope_data_records, only: lodiso
      implicit none (type, external)
      character (len = *), intent(in) :: isotope_file_path
      character(len=3), parameter  :: requested_status = "OLD"
      character(len=6), parameter  :: routine_name = "SETISO"
      character(len=9), parameter  :: record_format = "FORMATTED"
      character(len=14), parameter :: expected_file_header = "Atomic number:"
      integer           :: path_length, open_status, read_status
      logical           :: file_exists
      character(len=14) :: file_header
      inquire(file=isotope_file_path, exist=file_exists)
      if (.not.file_exists) then
         path_length = len_trim(isotope_file_path)
         write (stderr_unit,*)                                              &
                     routine_name,"- file: ",isotope_file_path(1:path_length)," does not exist"
         stop
      end if
      call openfl (22, isotope_file_path, record_format, requested_status, open_status)
      if (open_status /= 0) then
         write (stderr_unit, *) "Error opening isodata file: ", isotope_file_path(1:path_length)
         stop
      end if
      read (22, "(A)", iostat=read_status) file_header
      if (read_status/=0 .or. file_header/=expected_file_header) then
         write (stderr_unit, *) "Not an ISOtope Data File;"
         close(22)
         stop
      end if
      call lodiso
      close(22)

      return
      end subroutine setiso
end module isotope_data_input

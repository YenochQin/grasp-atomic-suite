! Read a case-insensitive yes/no response, prompting again for other characters.
module yes_no_input
    implicit none (type, external)
    private
    public :: getyn
contains
      logical function getyn ()
      use standard_io_units,     only: stdin_unit, stderr_unit
      implicit none (type, external)
      character (len = 1) :: response_character
    1 continue
      read (stdin_unit, "(A)") response_character
      if (response_character == "y" .or. response_character == "Y") then
         getyn = .true.
      else if (response_character == "n" .or. response_character == "N") then
         getyn = .false.
      else
         write (stderr_unit, *) "Expecting <y><cr> or <n><cr> ..."
         goto 1
      end if
      return
      end function getyn
end module yes_no_input

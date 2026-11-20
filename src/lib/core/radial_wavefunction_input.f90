! Check the radial wavefunction file signature on unit 23, load its orbital records, and close the file.
module radial_wavefunction_input
    implicit none (type, external)
    private
    public :: setrwfa
contains
      subroutine setrwfa(radial_file_path)
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use file_open, only: openfl
      use radial_wavefunction_records, only: lodrwf
      implicit none (type, external)
      character (len = *), intent(in) :: radial_file_path
      integer :: load_status, read_status
      character :: file_signature*6

      call openfl (23, radial_file_path, "UNFORMATTED", "OLD", load_status)
      if (load_status == 1) then
         write (stderr_unit, *) "Error when opening", radial_file_path(1:len_trim(radial_file_path))
         stop
      end if
      read (23, iostat=read_status) file_signature
      if (read_status/=0 .or. file_signature/="G92RWF") then
         write (stderr_unit, *) "This is not a Radial WaveFunction File;"
         close(23)
         stop
      end if
      call lodrwf (load_status)

      if (load_status /= 0) then
         write (stderr_unit, *) "Radial wavefunctions defined in CSL file", &
            " not found in Radial WaveFunction File"
         close(23)
         stop
      end if

      close(23)

      return
      end subroutine setrwfa
end module radial_wavefunction_input

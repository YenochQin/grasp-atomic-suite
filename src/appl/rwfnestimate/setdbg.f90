

!***********************************************************************
!                                                                      *
      subroutine setdbg(defnam)
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
!                                                                      *
!   Written by Farid A Parpia               Last update: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use calculation_input_defaults, only: inpci, ndef, ndump
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use file_open, only: openfl
      implicit none (type, external)
      character(len=*), intent(in) :: defnam
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ierr
      logical :: yes
      character :: filnam*256, form*11, status*3
!-----------------------------------------------
!
!
!   Initialise the arrays that control the debug printout
!
      ldbpa = .false.
!
      ldbpg = .false.
!
      ldbpr = .false.
!
      if (ndef == 0) return

      write (stderr_unit, *) "Generate debug printout?"
      yes = getyn()
      if (yes) then
!
!   The  .dbg  file is formatted; open it on unit 99
!
         form = "FORMATTED"
         status = "NEW"
!
         write (stderr_unit, *) "File ", defnam, " will be created as the ", &
            "ERWF DeBuG Printout File; "
         write (stderr_unit, *) "enter another file name if this is not ", &
            "acceptable; null otherwise:"
         read (*, "(A)") filnam
!
         if (len_trim(filnam) == 0) filnam = defnam
!
    4    continue
         call openfl (99, filnam, form, status, ierr)
         if (ierr /= 0) then
    5       continue
            write (stderr_unit, *) "Enter a name for the ERWF DeBuG Printout ", &
               "file that is to be created:"
            read (*, "(A)") filnam
            if (len_trim(filnam) == 0) goto 5
            goto 4
         end if
!
!   Set options for general printout
!
         write (stderr_unit, *) " Print out the machine constants used?"
         yes = getyn()
         if (yes) ldbpg(1) = .true.
         write (stderr_unit, *) " Print out the physical constants used?"
         yes = getyn()
         if (yes) ldbpg(2) = .true.
!
!   Set options for radial modules
!
         write (stderr_unit, *) " Printout from RADGRD?"
         yes = getyn()
         if (yes) ldbpr(1) = .true.
         write (stderr_unit, *) " Printout from NUCPOT?"
         yes = getyn()
         if (yes) ldbpr(2) = .true.
         write (stderr_unit, *) " Printout from TFPOT?"
         yes = getyn()
         if (yes) ldbpr(26) = .true.
!
!   Set options for angular modules
!
         write (stderr_unit, *) " Printout from LODCSL?"
         yes = getyn()
         if (yes) ldbpa(1) = .true.
!
      end if
!
      return
      end subroutine setdbg

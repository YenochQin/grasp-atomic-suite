!***********************************************************************
!                                                                      *
      subroutine setdbg (debug, file_path)
!                                                                      *
!   This routine opens the  .dbg  file and sets the arrays that con-   *
!   trol debug printout from the GENMCP program.                       *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
!                                                                      *
!   Written by Farid A Parpia               Last update: 08 Dec 1992   *
!   Modified by Xinghong He                 Last update: 03 Jul 1998   *
!
!   Used by mcpvu, mcpmpivu
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use calculation_input_defaults, only: ndef
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: debug
      character(len=*) :: file_path
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, open_status
      logical :: yes
!
      character(len = len (file_path)) :: filnam
      character(len=11) :: form
      character(len=3) :: status
!
!   Initialise the arrays that control the debug printout
!
      debug = .false.
      ibug1 = 0
      ibug2 = 0
      ibug3 = 0
      do i = 1, 5
         ldbpa(i) = .false.
         ldbpg(i) = .false.
      end do
!
      return        ! Skip these questions and return immediatel
!
      if (ndef == 0) then
         return
      end if

      write (stderr_unit,*) "Generate debug printout?"
      debug = getyn ()
      if (debug) then
!
!   The  .dbg  file is formatted; open it on unit 99
!
         form = "FORMATTED"
         status = "NEW"
!
         write (stderr_unit,*) "File ", file_path," will be created as the"    &
     &,                 " GENMCP DeBuG Printout File; "
         write (stderr_unit,*) " enter another file name if this is not "     &
     &,                 "acceptable; null otherwise:"
         read (*,"(A)") filnam
!
         if (len_trim (filnam) == 0) filnam = file_path
!
    4    call openfl (99, filnam, form, status, open_status)
         if (open_status /= 0) then
    5       write (stderr_unit,*) "Enter a name for the GENMCP DeBuG Printout"&
     &,                    " file that is to be created:"
            read (*,"(A)") filnam
            if (len_trim (filnam) == 0) goto 5
            goto 4
         end if
!
!   Set options for general printout
!
         write (stderr_unit,*) " Print out the machine constants used?"
         ldbpg(1) = getyn ()
!
!   Set options for angular modules
!
         write (stderr_unit,*) " Printout from LODCSL? (Not used)"
         ldbpa(1) = getyn ()
         write (stderr_unit,*) " Print out T coefficients?"
         ldbpa(2) = getyn ()
         write (stderr_unit,*) " Print out Coulomb V coefficients?"
         ldbpa(3) = getyn ()
         write (stderr_unit,*) " Print out sparse matrix definition arrays?"
         ldbpa(4) = getyn ()

      end if
!
      return
      end subroutine setdbg

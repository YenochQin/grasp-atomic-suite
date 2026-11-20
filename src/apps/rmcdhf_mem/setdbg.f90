!***********************************************************************
!                                                                      *
      subroutine setdbg(debug_file_path)
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 10 Dec 1992   *
!   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
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
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: debug_file_path
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: form = "FORMATTED"
      character(len=3), parameter :: status = "NEW"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, lendbg, lenfil, open_status
!GG      CHARACTER(LEN=120) :: FILNAM
      character(len = len (debug_file_path)) :: file_path
!-----------------------------------------------
!
!   Initialise the arrays that control the debug printout
!   These serve as the default settings.
!
      ldbpa = .false.
      ldbpg = .false.

      ldbpr = .false.
!
      if (ndef == 0) return

!   Even in non-default, the user can choose not to have debug
!   print-out

      write (stderr_unit, "(A)", advance="NO") "Generate debug output?  (y/n) "
      if (.not.getyn()) return

      lendbg = len_trim(debug_file_path)

      write (stderr_unit, *) "File  ", debug_file_path(1:lendbg), &
      "  will be                           created as the RSCF92 DeBuG Printout&
      & File;"
      write (stderr_unit, *) "enter another file name if this is not", &
         " acceptable; null otherwise:"

  123 continue
      read (*, "(A)") file_path

      file_path = adjustl(file_path)
      lenfil = len_trim(file_path)
      if (lenfil == 0) then
         file_path = debug_file_path
      else if (lenfil > lendbg) then
         write (stderr_unit, *) "File name too long, (> ", lendbg, "); redo..."
         goto 123
      end if

      call openfl (99, file_path, form, status, open_status)
      if (open_status /= 0) then
         write (stderr_unit, *) "File name not accepted; redo..."
         goto 123
      end if
!
!   Set options for general printout
!
      write (stderr_unit, *) "Print out the machine constants used?"
      ldbpg(1) = getyn()
      write (stderr_unit, *) "Print out the physical constants used?"
      ldbpg(2) = getyn()
      write (stderr_unit, *) "Printout from FNDBLK?"
      ldbpg(3) = getyn()
      write (stderr_unit, *) "Print out the Hamiltonian matrix?"
      ldbpg(4) = getyn()
      write (stderr_unit, *) "Print out the eigenvectors?"
      ldbpg(5) = getyn()
      ldbpg(1:5) = .true.
!
!   Set options for printout from radial modules
!
      write (stderr_unit, *) "Printout from RADGRD?"
      ldbpr(1) = getyn()
      write (stderr_unit, *) "Printout from NUCPOT?"
      ldbpr(2) = getyn()
      write (stderr_unit, *) "Printout from LODRWF?"
      ldbpr(3) = getyn()
      write (stderr_unit, *) "Print out I(ab) integrals?"
      ldbpr(4) = getyn()
      write (stderr_unit, *) "Print out Slater integrals?"
      ldbpr(10) = getyn()
      write (stderr_unit, *) "Make summary printout on progress", &
         " of each iteration in SOLVE?"
      ldbpr(22) = getyn()
      write (stderr_unit, *) "Tabulate and make printer plots", &
         " of subshell radial functions on", " each iteration in SOLVE?"
      ldbpr(23) = getyn()
      write (stderr_unit, *) "Tabulate and make printer plots", &
         " of subshell radial functions", " after each SCF cycle?"
      ldbpr(24) = getyn()
      write (stderr_unit, *) "Tabulate and make printer plots", &
         " of subshell radial functions on", " convergence?"
      ldbpr(25) = getyn()
      write (stderr_unit, *) "List compositions of exchange", " potentials?"
      ldbpr(27) = getyn()
      write (stderr_unit, *) "Tabulate and make printer plots", &
         " of exchange potentials?"
      ldbpr(28) = getyn()
      write (stderr_unit, *) "List compositions of direct", " potentials?"
      ldbpr(29) = getyn()
      write (stderr_unit, *) "Tabulate and make printer plots", &
         " of direct potentials?"
      ldbpr(30) = getyn()
      ldbpr(1:30) = .true.
!
!   Set options for printout of angular coefficients
!
      write (stderr_unit, *) " Printout from LODCSL?"
      ldbpa(1) = getyn()
      write (stderr_unit, *) " Print out T coefficients?"
      ldbpa(2) = getyn()
      write (stderr_unit, *) " Print out V coefficients?"
      ldbpa(3) = getyn()
      ldbpa(1:3) = .true.

      return
      end subroutine setdbg

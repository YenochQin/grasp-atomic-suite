!***********************************************************************
!                                                                      *
      subroutine setdbg(dbgfile)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:22:29   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use default_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use openfl_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: dbgfile
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character*9, parameter :: form = 'FORMATTED'
      character*3, parameter :: status = 'NEW'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, lendbg, lenfil, ierr
!GG      CHARACTER(LEN=120) :: FILNAM
      character(len = len (dbgfile)) :: filnam
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

      write (stderr_unit, '(A)') 'Generate debug output?  (y/n) '
      if (.not.getyn()) return

      lendbg = len_trim(dbgfile)

      write (stderr_unit, *) 'File  ', dbgfile(1:lendbg), &
      '  will be                           created as the RSCF92 DeBuG Printout&
      & File;'
      write (stderr_unit, *) 'enter another file name if this is not', &
         ' acceptable; null otherwise:'

  123 continue
      read (*, '(A)') filnam

      filnam = adjustl(filnam)
      lenfil = len_trim(filnam)
      if (lenfil == 0) then
         filnam = dbgfile
      else if (lenfil > lendbg) then
         write (stderr_unit,*) 'File name too long, (> ', lendbg, '); redo...'
         go to 123
      endif

      call openfl (99, filnam, form, status, ierr)
      if (ierr /= 0) then
         write (stderr_unit, *) 'File name not accepted; redo...'
         go to 123
      endif
!
!   Set options for general printout
!
      write (stderr_unit, *) 'Print out the machine constants used?'
      ldbpg(1) = getyn()
      write (stderr_unit, *) 'Print out the physical constants used?'
      ldbpg(2) = getyn()
      write (stderr_unit, *) 'Printout from FNDBLK?'
      ldbpg(3) = getyn()
      write (stderr_unit, *) 'Print out the Hamiltonian matrix?'
      ldbpg(4) = getyn()
      write (stderr_unit, *) 'Print out the eigenvectors?'
      ldbpg(5) = getyn()
!      LDBPG(1:5) = .TRUE.
!
!   Set options for printout from radial modules
!
      write (stderr_unit, *) 'Printout from RADGRD?'
      ldbpr(1) = getyn()
      write (stderr_unit, *) 'Printout from NUCPOT?'
      ldbpr(2) = getyn()
      write (stderr_unit, *) 'Printout from LODRWF?'
      ldbpr(3) = getyn()
      write (stderr_unit, *) 'Print out I(ab) integrals?'
      ldbpr(4) = getyn()
      write (stderr_unit, *) 'Print out Slater integrals?'
      ldbpr(10) = getyn()
      write (stderr_unit, *) 'Make summary printout on progress', &
         ' of each iteration in SOLVE?'
      ldbpr(22) = getyn()
      write (stderr_unit, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions on', ' each iteration in SOLVE?'
      ldbpr(23) = getyn()
      write (stderr_unit, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions', ' after each SCF cycle?'
      ldbpr(24) = getyn()
      write (stderr_unit, *) 'Tabulate and make printer plots', &
         ' of subshell radial functions on', ' convergence?'
      ldbpr(25) = getyn()
      write (stderr_unit, *) 'List compositions of exchange', ' potentials?'
      ldbpr(27) = getyn()
      write (stderr_unit, *) 'Tabulate and make printer plots', &
         ' of exchange potentials?'
      ldbpr(28) = getyn()
      write (stderr_unit, *) 'List compositions of direct', ' potentials?'
      ldbpr(29) = getyn()
      write (stderr_unit, *) 'Tabulate and make printer plots', &
         ' of direct potentials?'
      ldbpr(30) = getyn()
!      LDBPR(1:30) = .TRUE.
!
!   Set options for printout of angular coefficients
!
      write (stderr_unit, *) ' Printout from LODCSL?'
      ldbpa(1) = getyn()
      write (stderr_unit, *) ' Print out T coefficients?'
      ldbpa(2) = getyn()
      write (stderr_unit, *) ' Print out V coefficients?'
      ldbpa(3) = getyn()
!      LDBPA(1:3) = .TRUE.

      return
      end subroutine setdbg

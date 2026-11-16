!***********************************************************************
!                                                                      *
      subroutine setdbg (debug, fullname)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:11:16  12/23/06
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use default_c, only: ndef
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use openfl_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: debug
      character(len=*) :: fullname
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ierr
      logical :: yes
!
      character(len = len (fullname)) :: filnam
      character(len=11) :: form
      character(len=3) :: status
!
!   Initialise the arrays that control the debug printout
!
      do i = 1, 5
         ldbpa(i) = .false.
         ldbpg(i) = .false.
      enddo
!
      return      ! Skip these questions and return immediately
!
      if (ndef .eq. 0) then
         return
      endif

      write (stderr_unit,*) 'Generate debug printout?'
      debug = getyn ()
      if (debug) then
!
!   The  .dbg  file is formatted; open it on unit 99
!
         form = 'FORMATTED'
         status = 'NEW'
!
         write (stderr_unit,*) 'File ', fullname,' will be created as the'    &
     &,                 ' GENMCP DeBuG Printout File; '
         write (stderr_unit,*) ' enter another file name if this is not '     &
     &,                 'acceptable; null otherwise:'
         read (*,'(A)') filnam
!
         if (len_trim (filnam) .eq. 0) filnam = fullname
!
    4    call openfl (99, filnam, form, status, ierr)
         if (ierr .ne. 0) then
    5       write (stderr_unit,*) 'Enter a name for the GENMCP DeBuG Printout'&
     &,                    ' file that is to be created:'
            read (*,'(A)') filnam
            if (len_trim (filnam) .eq. 0) goto 5
            goto 4
         endif
!
!   Set options for general printout
!
         write (stderr_unit,*) ' Print out the machine constants used?'
         ldbpg(1) = getyn ()
!
!   Set options for angular modules
!
         write (stderr_unit,*) ' Printout from LODCSL? (Not used)'
         ldbpa(1) = getyn ()
         write (stderr_unit,*) ' Print out T coefficients?'
         ldbpa(2) = getyn ()
         write (stderr_unit,*) ' Print out Coulomb V coefficients?'
         ldbpa(3) = getyn ()
         write (stderr_unit,*) ' Print out sparse matrix definition arrays?'
         ldbpa(4) = getyn ()

      endif
!
      return
      end subroutine setdbg

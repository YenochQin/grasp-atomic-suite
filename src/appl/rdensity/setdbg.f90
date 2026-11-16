!***********************************************************************
!                                                                      *
      subroutine setdbg
!-----------------------------------------------
!                                                                      *
!   This subroutine sets the arrays that control debug printout from   *
!   the radial and angular modules of the GRASP92 suite.               *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
!                                                                      *
!   Written by Farid A Parpia               Last update: 24 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use debug_c
      use default_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use openfl_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: ierr
      logical   :: yes
      character :: filnam*256, defnam*11, form*11, status*3
!-----------------------------------------------
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

      write (6, *) 'Generate debug printout?'
      yes = getyn()
      if (yes) then
!
!   The  .dbg  file is formatted; open it on unit 99
!
         defnam = 'rdensity.dbg'
         form = 'FORMATTED'
         status = 'NEW'
!
         write (6, *) 'File  rdensity.dbg  will be created as the'
         write (6, *) ' rdensity DeBuG Printout File; enter another'
         write (6, *) ' file name if this is not acceptable;'
         write (6, *) ' null otherwise:'
         read (*, '(A)') filnam
!
         if (len_trim(filnam) == 0) filnam = defnam
!
    4    continue
         call openfl (99, filnam, form, status, ierr)
         if (ierr /= 0) then
    5       continue
            write (6, *) 'Enter a name for the rdensity DeBuG Printout'
            write (6, *) ' file that is to be created:'
            read (*, '(A)') filnam
            if (len_trim(filnam) == 0) go to 5
            go to 4
         endif
!
!   Set options for general printout
!
         write (6, *) ' Print out the machine constants used?'
         yes = getyn()
         if (yes) ldbpg(1) = .true.
         write (6, *) ' Print out the physical constants used?'
         yes = getyn()
         if (yes) ldbpg(2) = .true.
!
!   Set options for radial modules
!
         write (6, *) ' Printout from radial modules?'
         yes = getyn()
         if (yes) then
            write (6, *) ' Printout from RADGRD?'
            yes = getyn()
            if (yes) ldbpr(1) = .true.
            write (6, *) ' Printout from NUCPOT?'
            yes = getyn()
            if (yes) ldbpr(2) = .true.
            write (6, *) ' Printout from LODRWF?'
            yes = getyn()
            if (yes) ldbpr(3) = .true.
!
         endif
!
!   Set options for angular modules
!
         write (6, *) ' Printout from angular modules?'
         yes = getyn()
         if (yes) then
            write (6, *) ' Printout from LODCSL?'
            yes = getyn()
            if (yes) ldbpa(1) = .true.
            write (6, *) ' Print out V coefficients?'
            yes = getyn()
            if (yes) ldbpa(2) = .true.
         endif
!
      endif
!
      return
      end subroutine setdbg

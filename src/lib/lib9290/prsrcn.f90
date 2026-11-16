!***********************************************************************
!                                                                      *
      subroutine prsrcn(record, ncore, ioccs, ierr)
!                                                                      *
!   READs and parses a string that specifies a configuration.          *
!                                                                      *
!   Written by Farid A Parpia              Last revised: 16 Oct 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:15   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use orb_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: ncore
      integer, intent(out)          :: ierr
      character(len=256), intent(in) :: record
      integer, intent(out) :: ioccs(max_orbitals)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i,istart,iend,lenth,ios,npi,j,ishell,iostrt,ioend, &
         iocci,nkji
      character(len=2) :: symi
      character(len=5) :: form
      character(len=1), parameter :: cnum(3) = (/'1', '2', '3'/)
      character :: reci
!-----------------------------------------------
!
!   Initialise IOCCS for the peel subshells
!
      ioccs(ncore+1:orbital_count) = 0
!
!   Parse RECORD from left to right
!
      istart = 0
      i = 1
    2 continue
      reci = record(i:i)
      if (reci == '(') then
         iend = i - 1
    3    continue
         reci = record(iend:iend)
         if (reci == ' ') then
            iend = iend - 1
            go to 3
         else if (reci == '-') then
            read (record(iend-1:iend), '(1A2)') symi
            iend = iend - 2
         else
            symi(2:2) = ' '
            read (reci, '(1A1)') symi(1:1)
            iend = iend - 1
         endif
         lenth = iend - istart + 1
         form = '(1I'//cnum(lenth)//')'
         read (record(istart:iend), fmt=form, iostat=ios) npi
         if (ios /= 0) then
            write (6, *) 'PRSRCN: Principal quantum number ', record(istart:&
               iend)
            write (6, *) ' could not be decoded.'
            ierr = 1
            go to 6
         endif
         do j = ncore + 1, orbital_count
            if (orbital_principal_n(j)/=npi .or. orbital_symmetry_labels(j)/=symi) cycle
            ishell = j
            go to 5
         end do
         write (6, *) 'PRSRCL: Not a peel subshell.'
         ierr = 2
         go to 6
    5    continue
         iostrt = i + 1
      else if (reci == ')') then
         ioend = i - 1
         lenth = ioend - iostrt + 1
         form = '(1I'//cnum(lenth)//')'
         read (record(iostrt:ioend), fmt=form, iostat=ios) iocci
         if (ios /= 0) then
            write (6, *) 'PRSRCN: Occupation number ', record(iostrt:ioend)
            write (6, *) ' could not be decoded.'
            ierr = 3
            go to 6
         endif
         nkji = orbital_two_j(ishell)
         if (nkji <= 7) then
            if (iocci<0 .or. iocci>orbital_two_j(ishell)+1) then
               write (6, *) 'PRSRCN: Occupation specified'
               write (6, *) ' incorrectly for ', orbital_principal_n(ishell), orbital_symmetry_labels(ishell)
               write (6, *) ' subshell.'
               ierr = 4
               go to 6
            endif
         else
            if (iocci<0 .or. iocci>2) then
               write (6, *) 'PRSRCN: Occupation specified'
               write (6, *) ' incorrectly for ', orbital_principal_n(ishell), orbital_symmetry_labels(ishell)
               write (6, *) ' subshell.'
               ierr = 5
               go to 6
            endif
         endif
         ioccs(ishell) = iocci
         istart = 0
      else
         if (istart == 0) istart = i
      endif
!
      if (i < 256) then
         i = i + 1
         go to 2
      endif
!
      ierr = 0
!
    6 continue
      return
      end subroutine prsrcn

module getrsl_i
    implicit none (type, external)
    private
    public :: getrsl
contains
!***********************************************************************
!                                                                      *
      subroutine getrsl(indx, nsubs)
!                                                                      *
!   READs and parses a list of relativistic subshell labels delimit-   *
!   ed either by blanks or by commas. An asterisk may be used as the   *
!   `wildcard' character.                                              *
!
!   Output:
!     NSUBS - the # of orbitals parsed
!     indx(1:NSUBS) - indeces of these orbitals
!                                                                      *
!   Call(s) to: [LIB92]: LDIGIT.                                       *
!                                                                      *
!   Written by Farid A. Parpia             Last revised: 18 Dec 1992   *
!   Modified by Xinghong He                Last revised: 09 Jul 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orb_c
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ldigit_i, only: ldigit
      use convrt_i, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: nsubs
      integer, intent(inout) :: indx(*)
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istart, i, iend, nform, j, lenth, ios, npq, k
      logical :: found, nlany, nlok, npany, npok, nsany, nsok
      character :: nlq, nsq, reci, cnum*2, form*5, record*500
!-----------------------------------------------
!
!

      goto 2

    1 continue
      write (stderr_unit, *) " redo ..."
!
!   Read a record
!
    2 continue
      nsubs = 0
      read (*, "(A)") record

      write(734,"(a)") trim(record)
!
!   Parse RECORD from left to right
!
      istart = 0
      i = 1
    3 continue
      reci = record(i:i)
      if (reci/=" " .and. reci/=",") then
         if (istart == 0) istart = i
      else
         if (istart /= 0) then
            iend = i - 1
!
!   Parse the substring from left to right
!
!   (1) Determine the principal quantum number
!
            if (record(istart:istart) == "*") then
               npany = .true.
               istart = min(istart + 1,iend)
            else
               npany = .false.
               nform = 0
               do j = istart, iend
                  if (.not.ldigit(record(j:j))) cycle
                  nform = nform + 1
               end do
               if (nform == 0) then
                  write (stderr_unit, *) "GETRSL: Unable to interpret ", &
                     "the principal quantum number;"
                  goto 1
               end if
               call convrt (nform, cnum, lenth)
               form = "(1I"//cnum(1:lenth)//")"
               read (record(istart:istart+nform-1), form, iostat=ios) npq
               if (ios /= 0) then
                  write (stderr_unit, *) "GETRSL: Unable to interpret ", "string ", &
                     record(istart:iend-2), " as a principal quantum number"
                  goto 1
               end if
               istart = istart + nform
            end if
!
!   (2) Determine the orbital angular momentum quantum number
!
            nlq = record(istart:istart)
            if (nlq == "*") then
               nlany = .true.
            else
               nlany = .false.
            end if
!
!   (3) Determine the spin-orbit component
!
            if (iend > istart) then
               nsq = record(iend:iend)
               if (nsq == "*") then
                  nsany = .true.
               else if (nsq == "-") then
                  nsany = .false.
               else
                  write (stderr_unit, *) "GETRSL: Unable to interpret ", "string ", &
                     nsq, " as a spin-orbit component indicator"
                  goto 1
               end if
            else
               if (nlany) then
                  nsany = .true.
               else
                  nsany = .false.
                  nsq = " "
               end if
            end if
!
!
            found = .false.
            do j = 1, orbital_count
               npok = npany .or. orbital_principal_n(j)==npq
               nlok = nlany .or. nlq==orbital_symmetry_labels(j)(1:1)
               nsok = nsany .or. nsq==orbital_symmetry_labels(j)(2:2)
               if (.not.(npok .and. nlok .and. nsok)) cycle
               do k = 1, nsubs
                  if (indx(k) /= j) cycle
                  write (stderr_unit, *) "GETRSL: ", "Repeated subshell in list;"
                  goto 1
               end do
               found = .true.
               nsubs = nsubs + 1
               indx(nsubs) = j
            end do
!
            if (.not.found) then
               write (stderr_unit, *) "GETRSL: Subshell not occupied as ", &
                  " according to CSL  File;"
               goto 1
            end if
!
            istart = 0
         end if
      end if
!
      if (i < 500) then
         i = i + 1
         goto 3
      end if
!
      return
      end subroutine getrsl
end module getrsl_i

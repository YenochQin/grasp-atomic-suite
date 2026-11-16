module prsrsl_i
    implicit none (type, external)
    private
    public :: prsrsl
contains
!***********************************************************************
!                                                                      *
      subroutine prsrsl(nfile, id)
!                                                                      *
!   READs and parses a list of subshell labels on unit NFILE to load   *
!   COMMON  blocks  /ORB4/,  /ORB5/, and /ORB10/; the value of NW in   *
!   COMMON/ORB2/ is incremnted by 1 for each new subshell label; the   *
!   labels are delimited either by blanks or commas.                   *
!                                                                      *
!   Call(s) to: [LIB92] CONVRT.                                        *
!                                                                      *
!   Written by Farid A Parpia               Last revised 14 Sep 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_orbitals
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nfile
      integer , intent(in) :: id
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: istart, i, iend, ios, ii, iib2, lteger
      logical :: newrec
      character(len=1)   :: reci
      character(len=2)   :: cteger, sym
      character(len=6)   :: form
!      CHARACTER(LEN=500) :: RECORD, RECORD2
      character(len=1000) :: record, record2
      character(len=2) :: symlst(19)
!
      data symlst/ "s ", "p-", "p ", "d-", "d ", "f-", "f ", "g-", "g ", "h-", &
         "h ", "i-", "i ", "k-", "k ", "l-", "l ", "m-", "m"/
!-----------------------------------------------
      ios = 0
!
!   Read the records
!
      newrec = .false.
      read (nfile, "(A)") record
      if (id == 2) then
         read (nfile, "(A)") record2
         if (record2(1:3) == "CSF") then
            backspace (unit=nfile)
         else
            newrec = .true.
         end if
      end if
!
!   Parse RECORD from left to right
!
      istart = 0
      i = 1
    1 continue
      reci = record(i:i)
      if (reci/=" " .and. reci/=",") then
         if (istart == 0) istart = i
      else
         if (istart /= 0) then
            iend = i - 1
            reci = record(iend:iend)
            if (reci == "-") then
!              READ (RECORD(IEND-1:IEND),'(1A2)',IOSTAT = IOS) SYM
               sym = record(iend-1:iend)
               if (ios /= 0) then
                  write (stderr_unit, *) "PRSRSL: Symmetry ", record(iend-1:iend), &
                     " could not be decoded."
                  stop
               end if
               iend = iend - 2
            else
               sym(2:2) = " "
!              READ (RECI,'(1A1)',IOSTAT = IOS) SYM(1:1)
               sym(1:1) = reci
               if (ios /= 0) then
                  write (stderr_unit, *) "PRSRSL: Symmetry ", reci, " could not", &
                     " be decoded."
                  stop
               end if
               iend = iend - 1
            end if
            do ii = 1, 19
               if (sym /= symlst(ii)) cycle
               orbital_count = orbital_count + 1
               if (orbital_count > max_orbitals) then
                  write (stderr_unit, *) "PRSRSL: Number of subshells ", &
                     "exceeds allocation: plant NW was set to", max_orbitals
                  stop
               end if
               orbital_symmetry_labels(orbital_count) = sym
               iib2 = ii/2
               orbital_l(orbital_count) = iib2
               if (mod(ii,2) == 1) then
                  orbital_kappa(orbital_count) = (-iib2) - 1
                  orbital_two_j(orbital_count) = ii
               else
                  orbital_kappa(orbital_count) = iib2
                  orbital_two_j(orbital_count) = ii - 1
               end if
               call convrt (iend - istart + 1, cteger, lteger)
               form = "(1I"//cteger(1:lteger)//")"
               read (record(istart:iend), fmt=form, iostat=ios) orbital_principal_n(orbital_count)
               if (ios /= 0) then
                  write (stderr_unit, *) "PRSRSL: Principal quantum number ", record(&
                     istart:iend), " could not be decoded."
                  stop
               end if
               goto 3
            end do
            write (stderr_unit, *) "PRSRSL: Symmetry ", sym, " could not be decoded."
            stop
    3       continue
            istart = 0
         end if
      end if
!
!      IF (I < 500) THEN
      if (i < 1000) then
         i = i + 1
         goto 1
      end if
      if (newrec) then
         record = record2
         newrec = .false.
         i = 1
         istart = 0
         goto 1
      end if
!
      return
      end subroutine prsrsl
end module prsrsl_i

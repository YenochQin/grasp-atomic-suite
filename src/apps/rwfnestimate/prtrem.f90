

!***********************************************************************
!                                                                      *
      subroutine prtrem(all)
!                                                                      *
!   Prints a list of subshells that remain to be estimated.  ALL  is   *
!   .TRUE.  if all subshells have been estimated and  .FALSE. other-   *
!   wise.                                                              *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use initial_orbital_selection, only: max_orbitals, set
      use orbital_csf_basis, only: orbital_count,  orbital_principal_n, orbital_kappa, orbital_symmetry_labels
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical , intent(out) :: all
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /ORB2/
!     COMMON /ORB2/ NCF, NW, PNTRIQ
!     REAL   PNTRIQ
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: orbital_index, iend, ibeg, text_length
      character(len=12) :: number_text
      character(len=80) :: record*80
!-----------------------------------------------
!
!     Print *, ' Entering PRTREM'
!   Determine if there are any subshell radial wavefunctions that
!   remain to be estimated
!
      do orbital_index = 1, orbital_count
         if (set(orbital_index)) cycle
         all = .false.
         goto 2
      end do
      all = .true.
      goto 4
!
!   Print a list of subshell radial wavefunctions that remain to
!   be estimated; this list is no more than 80 characters wide
!
    2 continue
      write (stderr_unit, *) "The following subshell radial wavefunctions ", &
         "remain to be estimated:"
      iend = 0
      do orbital_index = 1, orbital_count
         if (set(orbital_index)) cycle
         if (iend > 75) then
            write (stderr_unit, *) record(1:iend)
            iend = 0
         end if
         if (iend > 0) then
            ibeg = iend + 1
            iend = ibeg
            record(ibeg:iend) = " "
         end if
         ibeg = iend + 1
         call convrt (orbital_principal_n(orbital_index), number_text, text_length)
         iend = ibeg + text_length - 1
         record(ibeg:iend) = number_text(1:text_length)
         if (orbital_kappa(orbital_index) < 0) then
            text_length = 1
         else
            text_length = 2
         end if
         ibeg = iend + 1
         iend = ibeg + text_length - 1
         record(ibeg:iend) = orbital_symmetry_labels(orbital_index)(1:text_length)
      end do
      if (iend > 1) write (stderr_unit, *) record(1:iend)
!
    4 continue
      return
      end subroutine prtrem

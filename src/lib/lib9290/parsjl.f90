!***********************************************************************
!                                                                      *
      subroutine parsjl(mode, ncore, record, loc, jx, njx, ierr)
!                                                                      *
!   READs and  parses a string that specifies angular momentum quan-   *
!   tum numbers.                                                       *
!                                                                      *
!   Call(s) to: [LIB92] CONVRT.                                        *
!                                                                      *
!   Written by Farid A Parpia              Last revised: 21 Dec 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:09   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i
      use orb_c, only: csf_count, orbital_count, csf_subshell_occupations
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: mode, ncore, loc
      integer, intent(out)          :: njx, ierr
      character(len=256), intent(in) :: record
      integer, intent(out) :: jx(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: njxmax, istart, i, iloop, ifrac, iend, lteger, j
      character(len=1) :: reci
      character(len=2) :: cteger
      character(len=6) :: form
!-----------------------------------------------
!
! This subroutine is here to replace the existing one which has been
! renamed as PARSJL_OLD. The purpose is to remove illegal GOTO into an
! IF block. The strategy is simple: copy the kernl to the bottom
! and add a do-loop over I.
! To restore (re-use) this subroutine, give the existing PARSJL a new
! name and then change the name of this subroutine back to PARSJL, and
! compile.
! XHH 1997.01.28
!
!
!   There cannot be more than JXMAX angular momenta specified; if
!   MODE is 1, the subshell quantum numbers are being read; if MODE
!   is 2, the intermediate and final angular momenta are being read
!
      if (mode == 1) then
         njxmax = 2*(orbital_count - ncore)
      else
         njxmax = orbital_count - ncore
      endif
!
!   Initialise NJX
!
      njx = 0
!
!   Parse RECORD from left to right
!
      istart = 0
      i = 1

! The original algorithm goes through the whole subroutine at least
! once, whatever the value of LOC. Thus we define another integer
! iloop to achieve this.
! XHH 1997.01.28

      iloop = max(1,loc)
      do i = 1, iloop
         reci = record(i:i)
         if (reci/=' ' .and. reci/=',' .and. reci/=';') then
            if (istart == 0) then
               istart = i
               ifrac = 0
            else
               if (reci == '/') then
                 ifrac = i
               endif
            endif
         else
            if (istart /= 0) then
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
               njx = njx + 1
               if (njx > njxmax) then
                  write (6, *) 'PARSJL: Too many angular momentum', &
                     ' quantum numbers specified;'
                  ierr = 1
                  go to 3
               endif
               iend = i - 1
               if (ifrac == 0) then
                  call convrt (iend - istart + 1, cteger, lteger)
                  form = '(1I'//cteger(1:lteger)//')'
                  read (record(istart:iend), fmt=form) j
                  if (j < 0) then
                     write (6, *) 'PARSJL: Negative angular momentum', &
                        ' quantum number found;'
                     ierr = 2
                     go to 3
                  endif
                  jx(njx) = 2*j
               else
                  call convrt (iend - ifrac, cteger, lteger)
                  form = '(1I'//cteger(1:lteger)//')'
                  read (record(ifrac+1:iend), fmt=form) j
                  if (j /= 2) then
                     write (6, *) 'PARSJL: The denominator of a', &
                        ' fractional quantum number must be 2;'
                     ierr = 3
                     go to 3
                  endif
                  call convrt (ifrac - istart, cteger, lteger)
                  form = '(1I'//cteger(1:lteger)//')'
                  read (record(istart:ifrac-1), fmt=form) j
                  if (j < 0) then
                     write (6, *) 'PARSJL: Negative angular momentum', &
                        ' quantum number found;'
                     ierr = 4
                     go to 3
                  endif
                  jx(njx) = j
               endif
               istart = 0
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            endif
         endif
      end do
!
      if (loc >= 1) then

! The following was accessed one extra time only when 1 <= I = LOC+1 .
! After the do-loop above, the value of I would be either LOC+1 or
! 2, depending on if LOC is greater than 1 or not, respectively


!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The following is exactly the same as those above
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         njx = njx + 1
         if (njx > njxmax) then
            write (6, *) 'PARSJL: Too many angular momentum', &
               ' quantum numbers specified;'
            ierr = 1
            go to 3
         endif
         iend = i - 1
         if (ifrac == 0) then
            call convrt (iend - istart + 1, cteger, lteger)
            form = '(1I'//cteger(1:lteger)//')'
            read (record(istart:iend), fmt=form) j
            if (j < 0) then
               write (6, *) 'PARSJL: Negative angular momentum', &
                  ' quantum number found;'
               ierr = 2
               go to 3
            endif
            jx(njx) = 2*j
         else
            call convrt (iend - ifrac, cteger, lteger)
            form = '(1I'//cteger(1:lteger)//')'
            read (record(ifrac+1:iend), fmt=form) j
            if (j /= 2) then
               write (6, *) 'PARSJL: The denominator of a', &
                  ' fractional quantum number must be 2;'
               ierr = 3
               go to 3
            endif
            call convrt (ifrac - istart, cteger, lteger)
            form = '(1I'//cteger(1:lteger)//')'
            read (record(istart:ifrac-1), fmt=form) j
            if (j < 0) then
               write (6, *) 'PARSJL: Negative angular momentum', &
                  ' quantum number found;'
               ierr = 4
               go to 3
            endif
            jx(njx) = j
         endif
         istart = 0
      endif
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      ierr = 0
!
    3 continue
      return
      end subroutine parsjl

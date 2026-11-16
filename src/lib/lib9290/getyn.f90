!***********************************************************************
!                                                                      *
      logical function getyn ()
!                                                                      *
!   This  subprogram reads a response on  the default input unit; the  *
!   response must be either 'y' or 'n'. GETYN is .TRUE. if 'y' is en-  *
!   tered and .FALSE. if 'n' is entered.                               *
!                                                                      *
!   Written by Farid A Parpia               Last update: 27 Aug 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:16   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c,     only: istdi, istde
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character (len = 1) :: rspns
!-----------------------------------------------
    1 continue
      read (istdi, '(A)') rspns
      if (rspns == 'y' .or. rspns == 'Y') then
         getyn = .true.
      else if (rspns == 'n' .or. rspns == 'N') then
         getyn = .false.
      else
         write (istde, *) 'Expecting <y><cr> or <n><cr> ...'
         go to 1
      endif
!
      return
      end function getyn

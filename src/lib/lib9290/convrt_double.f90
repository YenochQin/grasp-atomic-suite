!***********************************************************************
!                                                                      *
      subroutine convrt_double(intnum, cnum, lenth)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by G. Gaigalas,                                 May 2011   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E 10:46:53 2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!====================================================================

!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,   intent(in)    :: intnum
      integer,   intent(out)   :: lenth
      character, intent(inout) :: cnum*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                           :: intnumgg
      character(len=6)                  :: form
      character(len=2) :: c1020(0:10)
      character           :: c19(9)
!
      data c19 /'1','2','3','4','5','6','7','8','9'/
      data c1020 /'10','11','12','13','14','15','16','17','18','19','20'/
!-----------------------------------------------
      if(mod(intnum,2) == 0) then
         intnumgg = intnum/2
      else
         intnumgg = intnum
      endif
!
      if (intnumgg < 0) then
         lenth = log10(dble((-intnumgg))) + 2
      else if (intnumgg == 0) then
         lenth = 1
      else
         lenth = log10(dble(intnumgg)) + 1
      endif
!
!   Ensure that the length of CNUM as dimensioned is adequate;
!   stop with an error message if it isn't
!
      if (lenth > len(cnum)) then
         write (6, *) 'CONVRT_DOUBLE: Length of CNUM inadeuate.'
         error stop
      else
         if (lenth <= 9) then
            form = '(1I'//c19(lenth)//')'
            write (cnum(1:lenth), form(1:5)) intnumgg
         else
            form = '(1I'//c1020(lenth-10)//')'
            write (cnum(1:lenth), form(1:6)) intnumgg
         endif
         if(mod(intnum,2) /= 0) then
            if (lenth+2 > len(cnum)) then
               write (6, *) 'CONVRT_DOUBLE: Length of CNUM inadeuate.'
               error stop
            else
               cnum(1:lenth+2) = cnum(1:lenth)//'/2'
               lenth = lenth + 2
            endif
         endif
      endif
      return
      end subroutine convrt_double

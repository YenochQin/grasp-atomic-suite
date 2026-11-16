module convrt_i
    implicit none (type, external)
    private
    public :: convrt
contains
!***********************************************************************
!                                                                      *
      subroutine convrt(intnum, cnum, lenth)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 22 Sep 1992   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
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
      character(len=6)                  :: form
      character(len=2) :: c1020(0:10)
      character           :: c19(9)
!
      data c19 /"1","2","3","4","5","6","7","8","9"/
      data c1020 /"10","11","12","13","14","15","16","17","18","19","20"/
!-----------------------------------------------
      if (intnum < 0) then
         lenth = log10(dble((-intnum))) + 2
      else if (intnum == 0) then
         lenth = 1
      else
         lenth = log10(dble(intnum)) + 1
      end if
!
!   Ensure that the length of CNUM as dimensioned is adequate;
!   stop with an error message if it isn't
!
      if (lenth > len(cnum)) then
         write (6, *) "CONVRT: Length of CNUM inadeuate."
         error stop
      else
         if (lenth <= 9) then
            form = "(1I"//c19(lenth)//")"
            write (cnum(1:lenth), form(1:5)) intnum
         else
            form = "(1I"//c1020(lenth-10)//")"
            write (cnum(1:lenth), form(1:6)) intnum
         end if
      end if
!
      return
      end subroutine convrt
end module convrt_i

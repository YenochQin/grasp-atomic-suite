module convrt2_i
    implicit none (type, external)
    private
    public :: convrt2
contains
!***********************************************************************
!                                                                      *
      subroutine convrt2(intnum, cnum, lenth, from)
!                                                                      *
!   Converts the  INTEGER number  INTNUM  into the  CHARACTER string   *
!   CNUM of length LENTH. INTEGER lengths of up to 64 bits are acco-   *
!   modated.                                                           *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 22 Sep 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: intnum
      integer, intent(out) :: lenth
      character, intent(inout) :: cnum*(*)
      character, intent(in) :: from*(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character :: form*6
      character :: c1020(0:10)*2
      character :: c19(9)
!
!
      data c19/ "1", "2", "3", "4", "5", "6", "7", "8", "9"/
      data c1020/ "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", &
         "20"/
!
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
         write (stderr_unit, *) "CONVRT: Length of CNUM inadeuate. (from:", from, ")"
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
      end subroutine convrt2
end module convrt2_i

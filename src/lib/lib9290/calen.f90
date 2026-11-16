module calen_i
    implicit none (type, external)
    private
    public :: calen
contains
!***********************************************************************
!                                                                      *
      subroutine calen(jtime, jdate)
!                                                                      *
!   Loads the character strings JTIME and JDATE with the time of day   *
!   and the date when called.                                          *
!                                                                      *
!                                         Last revision: 25 Sep 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=10), intent(out) :: jtime
      character(len=8), intent(out) :: jdate

      call date_and_time(jdate, jtime)

      return
      end subroutine calen
end module calen_i

module ldigit_i
    implicit none (type, external)
    private
    public :: ldigit
contains
!***********************************************************************
!                                                                      *
      logical function ldigit (cst)
!                                                                      *
!   .TRUE.  if  CST  is the ASCII representation of a decimal digit;   *
!   .FALSE. otherwise.                                                 *
!                                                                      *
!   Written by Farid A. Parpia             Last revised: 16 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: cst
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      character :: cdgt(0:9)
!-----------------------------------------------
!
!
!
      data cdgt/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
!
      do i = 0, 9
         if (cst /= cdgt(i)) cycle
         ldigit = .true.
         go to 2
      end do
      ldigit = .false.
!
    2 continue
      return
      end function ldigit
end module ldigit_i

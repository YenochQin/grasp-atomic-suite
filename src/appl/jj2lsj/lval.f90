!
!***********************************************************************
!                                                                      *
      integer function lval (symbol)
!                                                                      *
!   Modified by Gediminas Gaigalas,                   September 1997   *
!   Modified by G. Gaigalas,                                May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01
!...Switches:
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: symbol
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer           :: locate
      character(len=26) :: set
!
      data set/ 'spdfghiklmnoqSPDFGHIKLMNOQ'/
!-----------------------------------------------
      locate = index(set,symbol)
      if (locate <= 13) then
         lval = locate - 1
      else
         lval = locate - 14
      endif
      return
      end function lval

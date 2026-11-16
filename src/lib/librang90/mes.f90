!*******************************************************************
!                                                                  *
      subroutine mes(i)
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Vilnius,  Lithuania                             December 1993  *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      implicit none (type, external)
      integer, intent(in) :: i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j
      character(len=10) :: string5(6)
!
      data string5/'   I T L S',' I T L S 2',' I T L S 3', &
                   'AW P 1 L S','WA P 1 L S','       W 1'/
!-----------------------------------------------
      if(i > 50) then
        j=i-50
        write(6,'(A)') ' error in func./sub. '
        write(6,'(20X,A10)') string5(j)
        write(6,'(A)') ' susimaise f sluoksnio termu kodavimas  '
      else
        write(6,'(A)') ' yra daugiau nei 2 ele. sluoks. f,g,h,i,k,l,m'
        write(6,'(3X,I5)') i
        if(i == 1) then
          write(6,'(A)') ' error in Subroutine   W 1 G '
        elseif(i == 2) then
          write(6,'(A)') ' error in Subroutine   W 1 '
        elseif(i == 11) then
          write(6,'(A)') ' error in Function     I T L S  '
        elseif(i == 12) then
          write(6,'(A)') ' error in Function     I T L S 2  '
        elseif(i == 13) then
          write(6,'(A)') ' error in Function     I T L S 3  '
        elseif(i == 30) then
          write(6,'(A)') ' error in Subroutine   A 1 A 2 A 3 A 4 L S '
        elseif(i == 31) then
          write(6,'(A)') ' error in Subroutine   A 1 A 2 L S '
        elseif(i == 32) then
          write(6,'(A)') ' error in Subroutine   A 1 A 2 W 3 L S '
        elseif(i == 33) then
          write(6,'(A)') ' error in Subroutine   A 1 A W 2 L S '
        elseif(i == 34) then
          write(6,'(A)') ' error in Subroutine   W A 1 A 2 L S '
        elseif(i == 35) then
          write(6,'(A)') ' error in Subroutine   W 1 W 2 L S '
        else
          write(6,'(A)') ' error in unknown Subroutine  '
        endif
      endif
      write(6,'(A)') ' Contact to   G. Gediminas please ! ! ! '
      stop
      end subroutine mes

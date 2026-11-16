!*******************************************************************
!                                                                  *
      subroutine el4(jja,jjb,ja,jb,jc,jd,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!     SUBROUTINE CALLED: EL41                                      *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use m_c,            only: npeel
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use el41_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,jc,jd,icolbrei
!-----------------------------------------------
      if(npeel <= 2)return
      if(ja == jb) then
        call el41(jja,jjb,jc,jd,ja,1,ja,jb,jc,jd,icolbrei)
      else if(jc == jd) then
        call el41(jja,jjb,ja,jb,jc,2,ja,jb,jc,jd,icolbrei)
      else
        write(99,100)
        stop
      end if
      return
  100 format(5x,'ERRO IN EL4 ')
      end subroutine el4

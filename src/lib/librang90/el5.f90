module el5_i
    implicit none (type, external)
    private
    public :: el5
contains
!*******************************************************************
!                                                                  *
      subroutine el5(jja,jjb,ja,jb,jc,jd,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!      SUBROUTINE CALLED: EL51,EL52,EL53                           *
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
      use m_c,          only: npeel
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use el51_i, only: el51
      use el52_i, only: el52
      use el53_i, only: el53
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,jc,jd,icolbrei
!-----------------------------------------------
      if(npeel <= 3)return
      if(jb < jc) then
        call el51(jja,jjb,ja,jb,jc,jd,1,icolbrei)
      else if(ja > jd.and.jb > jd) then
        call el51(jja,jjb,jc,jd,ja,jb,2,icolbrei)
      else if(jb > jc.and.jb < jd.and.ja < jc) then
        call el52(jja,jjb,ja,jc,jb,jd,1,icolbrei)
      else if(jb > jc.and.jb > jd.and.ja > jc) then
        call el52(jja,jjb,jc,ja,jd,jb,2,icolbrei)
      else if(jb > jc.and.jb > jd.and.ja < jc) then
        call el53(jja,jjb,ja,jc,jd,jb,1,icolbrei)
      else if(jb > jc.and.jb < jd.and.ja > jc) then
        call el53(jja,jjb,jc,ja,jb,jd,2,icolbrei)
      else
        write(99,100)
        stop
      end if
      return
  100 format(5x,"ERRO IN EL5 ")
      end subroutine el5
end module el5_i

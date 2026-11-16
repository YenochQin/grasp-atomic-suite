!*******************************************************************
!                                                                  *
      subroutine eile(ja,jb,jc,jaa,jbb,jcc)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 02  ------------   *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ja, jb, jc
      integer, intent(out) :: jaa, jbb, jcc
!-----------------------------------------------
      jaa=ja
      jcc=ja
      if(jaa > jb)jaa=jb
      if(jcc < jb)jcc=jb
      if(jaa > jc)jaa=jc
      if(jcc < jc)jcc=jc
      if((ja > jaa).and.(ja < jcc))jbb=ja
      if((jb > jaa).and.(jb < jcc))jbb=jb
      if((jc > jaa).and.(jc < jcc))jbb=jc
      return
      end subroutine eile

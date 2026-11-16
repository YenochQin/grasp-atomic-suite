module el3_i
    implicit none (type, external)
    private
    public :: el3
contains
!*******************************************************************
!                                                                  *
      subroutine el3(jja,jjb,ja,jb,jc,jd,icolbrei)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     SUBROUTINE CALLED: EL31,EL32,EL33                            *
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
      use m_c,         only:  npeel
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use el31_i
      use el32_i
      use el33_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jja,jjb,ja,jb,jc,jd,icolbrei
!-----------------------------------------------
      if(npeel <= 1)return
      if(jb == jd) then
        if(ja == jb.or.jc == jb) then
          if(ja == jc)go to 10
          if(jc /= jb) then
            call el32(jja,jjb,jc,ja,ja,jb,jc,jd,icolbrei)
          else
            call el31(jja,jjb,jc,ja,ja,jb,jc,jd,icolbrei)
          end if
        else
          call el33(jja,jjb,jc,ja,jb,1,ja,jb,jc,jd,icolbrei)
        end if
        return
      else if(ja == jc) then
        if(jb == ja.or.jd == ja) then
          if(jb == jd)go to 10
          if(jd /= ja) then
            call el32(jja,jjb,jd,jb,ja,jb,jc,jd,icolbrei)
          else
            call el31(jja,jjb,jd,jb,ja,jb,jc,jd,icolbrei)
          end if
        else
          call el33(jja,jjb,jd,jb,ja,1,ja,jb,jc,jd,icolbrei)
        end if
        return
      else if(ja == jd) then
        if(jb == ja.or.jc == ja) then
          if(jb == jc)go to 10
          if(jc /= jd) then
             call el32(jja,jjb,jc,jb,ja,jb,jc,jd,icolbrei)
          else
             call el31(jja,jjb,jc,jb,ja,jb,jc,jd,icolbrei)
          end if
        else
          call el33(jja,jjb,jc,jb,ja,2,ja,jb,jd,jc,icolbrei)
        end if
        return
      else if(jb == jc) then
        if(ja == jb.or.jd == jb) then
          if(ja == jd)go to 10
          if(jd /= jb) then
            call el32(jja,jjb,jd,ja,ja,jb,jc,jd,icolbrei)
          else
            call el31(jja,jjb,jd,ja,ja,jb,jc,jd,icolbrei)
          end if
        else
          call el33(jja,jjb,jd,ja,jb,2,ja,jb,jd,jc,icolbrei)
        end if
        return
      end if
   10 write(99,100)
  100 format(5x,'ERRO IN EL3  PMGG RAGG')
      stop
      end subroutine el3
end module el3_i

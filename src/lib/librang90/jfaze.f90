module jfaze_i
    implicit none (type, external)
    private
    public :: jfaze
contains
!*******************************************************************
!                                                                  *
      integer function jfaze(i1,i2,i3,i4)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 20  -------------  *
!                                                                  *
!     DETERMINATE THE PHASE FACTOR WHICH APPEAR FROM PERMUTATION   *
!     OF OPERATORS OF SECOND QUANTIZATION                          *
!                                                                  *
!     NO FUNCTION CALLED                                           *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i1, i2, i3, i4
!-----------------------------------------------
      jfaze=1
      if(i1 > i2)jfaze=-jfaze
      if(i1 > i3)jfaze=-jfaze
      if(i1 > i4)jfaze=-jfaze
      if(i2 > i3)jfaze=-jfaze
      if(i2 > i4)jfaze=-jfaze
      if(i3 > i4)jfaze=-jfaze
      return
      end function jfaze
end module jfaze_i

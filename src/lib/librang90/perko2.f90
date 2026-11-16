!*******************************************************************
!                                                                  *
      subroutine perko2(ja1,ja2,ja3,ja4,i)
!                                                                  *
!   --------------  SECTION METWO    SUBPROGRAM 23  -------------  *
!                                                                  *
!     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
!                                               (GENERAL CASE)     *
!                                                                  *
!     SUBROUTINE CALLED: PERKO1                                    *
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
      use trk_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use perko1_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja1, ja2, ja3, ja4, i
!-----------------------------------------------
      call perko1(ja1,bk1,ik1,bd1,id1)
      if(i == 1)return
      call perko1(ja2,bk2,ik2,bd2,id2)
      if(i == 2)return
      call perko1(ja3,bk3,ik3,bd3,id3)
      if(i == 3)return
      call perko1(ja4,bk4,ik4,bd4,id4)
      return
      end subroutine perko2

module rmeajj9_i
    implicit none (type, external)
    private
    public :: rmeajj9
contains
!*******************************************************************
!                                                                  *
      subroutine rmeajj9(it,lq,j,its,lqs,j1s,coef)
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
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: it, lq, j, its, lqs, j1s
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jt, jts
      integer :: ids3(20),ids4(20),ids5(20),ids6(20),ids7(20),idv7(20),ids8(20),&
      idv8(20),ids9(20),ids10(20),ids11(20),idv11(20),ids12(20),idv12(20),ids13(20),idv13(20),ids14(20),   &
      ids15(20),idv15(20),ids16(20),idv16(20),ids17(20),idv17(20),ids18(20),idv18(20)
!-----------------------------------------------
      data ids3/2*0,660,0,1664,0,1650,130,0,816,0,1680,8*0/
      data ids4/2*0,-55770,-20592,6292,-13260,0,12740,-19110,     &
      79968,-47880,-25410,-19278,7*0/
      data ids5/0,-9009,-4056,-260,585,16575,0,25200,4200,-7140,  &
      -475,504,-14280,13566,-4250,5*0/
      data ids6/0,-4992,2028,0,-1920,0,-12870,7350,0,8160,        &
      0,1512,0,-3648,0,-9000,4*0/
      data ids7/-2184,-63,-59904,-302460,5265,848691,0,-145152,   &
      217728,1049580,287337,5184,261120,691866,-750,0,-76608,3*0/
      data idv7/253,23,19481,19481,161,253253,1,36179,            &
      36179,2*36179,299,3*36179,1,3289,3*1/
      data ids8/1224,2652,188598,-31824,204,-12996,0,25500,-38250,&
      -768,81396,3213,-60543,-3066144,727776,207,-41553,3*0/
      data idv8/1265,115,13915,2783,23,13915,1,3*2783,13915,115,  &
      2783,236555,47311,17,21505,3*1/
      data ids9/3380,-2340,5460,-1400,3150,-3570,0,9000,1500,2550,&
      -5320,-10710,1050,1140,6300,8550,-330,5750,2*0/
      data ids10/0,2160,6240,0,-32,0,-21450,-9610,0,2688,         &
      0,7140,0,20160,0,6156,2*0,-10164,0/
      data ids11/0,132,-52728,196,-50,-24990,0,-21160,31740,26250,&
      12920,-357,9583,-344988,5700,171,-15,-4830,-84,0/
      data idv11/1,5,6655,1331,11,1331,1,4*1331,55,1331,6655,     &
      1331,11,2*121,11,1/
      data ids12/2*0,-209950,77520,12920,-25688,0,-4522,2261,     &
      -48640,-285144,931,273885,-112908,-2138580,137781,6654375,  &
      -59616,284089,0/
      data idv12/2*1,2*9317,231,9317,1,3993,1331,22627,9317,44,   &
      90508,429913,158389,3740,156332,39083,17765,1/
      data ids13/2*0,1530,13056,29376,720,0,-1890,-315,101124,    &
      -4560,-13965,-35131,13500,-685900,-1197,-28875,-5060,759,0/
      data idv13/2*1,77,3*1001,1,2*143,2431,1001,572,9724,2431,   &
      17017,2*884,221,17,1/
      data ids14/4*0,22848,0,-121550,4590,0,-32832,0,45220,0,     &
      -31680,0,82764,2*0,144716,0/
      data ids15/5*0,2128,0,-1938,2907,-4860,-17136,-1309,-8505,  &
      15876,420,1287,-6075,-132,-253,-650/
      data idv15/5*1,143,1,3*143,2717,52,572,247,13,20,988,19,95, &
      19/
      data ids16/7*0,570,95,4104,504,-1463,-39501,-60516,-1596,   &
      3933,621,-840,-805,390/
      data idv16/7*1,2*13,221,65,260,884,1105,221,68,3740,2*17,   &
      11/
      data ids17/9*0,-5796,17664,-5313,1771,-16632,-88,693,       &
      -94269,192500,30030,24570/
      data idv17/9*1,221,1235,65,221,1615,2*17,1615,7429,323,437/
      data ids18/13*0,-15000,280,-1170,-48750,-15600,3510,-33930/
      data idv18/13*1,323,2*17,3553,437,19,253/
!
      coef=zero
      if(it < its) then
        jt=it-25
        jts=its-45
      else
        jt=its-25
        jts=it-45
      endif
      if(jts == 1) then
        if(jt == 7) coef=-dsqrt(dble(60))
      elseif(jts == 2) then
        if(jt == 8) coef=-dsqrt(dble(12))
        if(jt == 9) coef=-dsqrt(dble(8))
      elseif(jts == 3)then
        coef=-dsqrt(dble(ids3(jt))/dble(33))
      elseif(jts == 4)then
        if(ids4(jt) >= 0) then
           coef=dsqrt(dble(ids4(jt))/dble(3003))
        else
           coef=-dsqrt(-dble(ids4(jt))/dble(3003))
        endif
      elseif(jts == 5)then
        if(ids5(jt) >= 0) then
           coef=dsqrt(dble(ids5(jt))/dble(715))
        else
           coef=-dsqrt(-dble(ids5(jt))/dble(715))
        endif
      elseif(jts == 6)then
        if(ids6(jt) >= 0) then
           coef=dsqrt(dble(ids6(jt))/dble(143))
        else
           coef=-dsqrt(-dble(ids6(jt))/dble(143))
        endif
      elseif(jts == 7)then
        if(ids7(jt) >= 0) then
           coef=dsqrt(dble(ids7(jt))/dble(idv7(jt)))
        else
           coef=-dsqrt(-dble(ids7(jt))/dble(idv7(jt)))
        endif
      elseif(jts == 8)then
        if(ids8(jt) >= 0) then
           coef=dsqrt(dble(ids8(jt))/dble(idv8(jt)))
        else
           coef=-dsqrt(-dble(ids8(jt))/dble(idv8(jt)))
        endif
      elseif(jts == 9)then
        if(ids9(jt) >= 0) then
           coef=dsqrt(dble(ids9(jt))/dble(325))
        else
           coef=-dsqrt(-dble(ids9(jt))/dble(325))
        endif
      elseif(jts == 10)then
        if(ids10(jt) >= 0) then
           coef=dsqrt(dble(ids10(jt))/dble(165))
        else
           coef=-dsqrt(-dble(ids10(jt))/dble(165))
        endif
      elseif(jts == 11)then
        if(ids11(jt) >= 0) then
           coef=dsqrt(dble(ids11(jt))/dble(idv11(jt)))
        else
           coef=-dsqrt(-dble(ids11(jt))/dble(idv11(jt)))
        endif
      elseif(jts == 12)then
        if(ids12(jt) >= 0) then
           coef=dsqrt(dble(ids12(jt))/dble(idv12(jt)))
        else
           coef=-dsqrt(-dble(ids12(jt))/dble(idv12(jt)))
        endif
      elseif(jts == 13)then
        if(ids13(jt) >= 0) then
           coef=dsqrt(dble(ids13(jt))/dble(idv13(jt)))
        else
           coef=-dsqrt(-dble(ids13(jt))/dble(idv13(jt)))
        endif
      elseif(jts == 14)then
        if(ids14(jt) >= 0) then
           coef=dsqrt(dble(ids14(jt))/dble(715))
        else
           coef=-dsqrt(-dble(ids14(jt))/dble(715))
        endif
      elseif(jts == 15)then
        if(ids15(jt) >= 0) then
           coef=dsqrt(dble(ids15(jt))/dble(idv15(jt)))
        else
           coef=-dsqrt(-dble(ids15(jt))/dble(idv15(jt)))
        endif
      elseif(jts == 16)then
        if(ids16(jt) >= 0) then
           coef=dsqrt(dble(ids16(jt))/dble(idv16(jt)))
        else
           coef=-dsqrt(-dble(ids16(jt))/dble(idv16(jt)))
        endif
      elseif(jts == 17)then
        if(ids17(jt) >= 0) then
           coef=dsqrt(dble(ids17(jt))/dble(idv17(jt)))
        else
           coef=-dsqrt(-dble(ids17(jt))/dble(idv17(jt)))
        endif
      elseif(jts == 18)then
        if(ids18(jt) >= 0) then
           coef=dsqrt(dble(ids18(jt))/dble(idv18(jt)))
        else
           coef=-dsqrt(-dble(ids18(jt))/dble(idv18(jt)))
        endif
      else
        write(0,'(A,4I5)') ' IT ITS JT JTS= ',it,its,jt,jts
        write(0,'(A)') ' ERROR IN SUB. RMEAJJ9 '
        stop
      endif
      return
      end subroutine rmeajj9
end module rmeajj9_i

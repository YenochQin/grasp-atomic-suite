module rmew7bjj_i
    implicit none (type, external)
    private
    public :: rmew7bjj
contains
!*******************************************************************
!                                                                  *
      subroutine rmew7bjj(j1,j2,k1,k2,coef)
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
      use ribojj_c, only: imgnjj, imgtjj, impnjj, imptjj
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rumtjj_i, only: rumtjj
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: j1, j2, k1, k2
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ll, iq1, lv1, il1, iq2, lv2, il2, ji1, ji2, ifaz, j, l
      integer  :: ipr1(6)
      integer  :: ipr2(8)
      integer :: i12ps(21),i12pv(21),i03ps(21),i03pv(21),i14ps(21), &
      i14pv(21),i05ps(21),i05pv(21),i16ps(21),i16pv(21),i07ps(21),i07pv(21)
      integer :: i12as(36),i12av(36),i03as(36),i03av(36),i14as(36), &
      i14av(36),i05as(36),i05av(36),i16as(36),i16av(36),i07as(36),i07av(36)
!-----------------------------------------------
      data ipr1/0,5,9,12,14,15/
      data ipr2/0,7,13,18,22,25,27,28/
      data i12ps/-252,1056,144,3*0,507,-88,-325,2*0,200,-520,80,  &
      0,-6125,-560,0,390,-3840,2040/
      data i12pv/10,70,7,3*1,10,1,14,2*1,3,21,2*1,462,11,1,539,   &
      2*49/
      data i12as/0,-50,6*0,-1280,990,2640,1950,4*0,-480,4*0,-360, &
      1872,42,390,3*0,48,2*0,234,0,1040,340,0/
      data i12av/8*1,4*49,4*1,49,4*1,539,49,1,11,6*1,7,1,11,7,1/
      data i03ps/-1188,-196,0,-234,2*0,189,0,-1911,1470,0,-56,3*0,&
      394805,5250,53760,-78,17408,12920/
      data i03pv/70,10,1,7,2*1,110,1,242,121,5*1,22022,121,1573,  &
      2*847,1001/
      data i03as/8*0,110,0,-240,4*0,-1920,0,52,-224,2*0,32490,2*0,&
      -7644,0,12,224,2*0,-52,0,-2040,-364,0,1292/
      data i03av/8*1,7,1,7,4*1,77,1,7,11,2*1,847,2*1,121,1,70,10, &
      2*1,70,1,77,121,1,77/
      data i14ps/0,54,-528,546,-378,0,-4335,-96,11271,3822,0,120, &
      12000,-624,-960,-30345,-210,-228480,580476,146880,-627912/
      data i14pv/1,2*7,2*11,1,110,11,1694,121,2*1,77,2*11,3146,   &
      121,1573,2*5929,7007/
      data i14as/3*0,-90,4*0,2640,480,-360,20592,-42,390,2*0,     &
      468180,2*0,21840,0,-359424,-6750,0,917280,-10710,2*0,36,2*0,&
      -858,0,-5304,-69768,0/
      data i14av/8*1,2*49,2*539,1,11,2*1,5929,2*1,847,1,5929,539, &
      1,5929,121,2*1,11,2*1,7,1,121,847,1/
      data i05ps/3*0,144,14,0,975,0,-50421,336,-7000,-88,3*0,     &
      -6845,-3360,14280,-1836,-103360,-28424/
      data i05pv/3*1,7,2*1,22,1,2002,11,143,4*1,26026,143,1859,77,&
      1001,143143/
      data i05as/10*0,390,2*0,-70,3*0,32,14,2*0,-576,2*0,-210,0,  &
      63888,-154,0,-8568,-176,0,1938,1088,0,-28424/
      data i05av/10*1,7,6*1,7,3*1,77,2*1,11,1,1183,13,1,169,7,1,  &
      91,11,1,1183/
      data i16ps/3*0,576,390,-144,0,520,-49,-12480,-408,520,      &
      -1960,-1664,3264,552250,-43520,-38760,-2652,15504,38760/
      data i16pv/3*1,11,77,7,1,11,2*121,11,3,33,2*11,4719,847,    &
      11011,2*121,143/
      data i16as/6*0,-130,3*0,390,-48,-234,1040,340,0,3120,2*0,   &
      -1170,0,18720,-36,858,-5304,-69768,2*0,1020,4*0,-33592,     &
      31654,0/
      data i16av/10*1,11,1,7,11,7,1,121,2*1,121,1,121,11,7,121,   &
      847,2*1,11,4*1,2*121,1/
      data i07ps/4*0,162,272,2*0,11025,4624,-1632,-120,3*0,       &
      1224510,306000,12558240,-6460,77520,-297160/
      data i07pv/4*1,2*7,2*1,1573,121,143,4*1,20449,11011,143143, &
      121,1573,1859/
      data i07as/13*0,60,4*0,1600,0,2040,4410,2*0,18360,0,-18816, &
      -11016,0,11628,34,0,7752,9690,0,222870/
      data i07av/18*1,77,1,77,121,2*1,121,1,845,455,1,1183,5,1,   &
      143,121,1,1859/
!
      coef = zero
      ll = 7
      call rumtjj(j1,ll,iq1,lv1,il1)
      call rumtjj(j2,ll,iq2,lv2,il2)
      if(j1 > j2) then
        ji1=j2
        ji2=j1
        ifaz=il2-il1+iq2-iq1
        else
          ji1=j1
          ji2=j2
          ifaz=4
      end if
      if(j1 > 17) then
        ji1=ji1-17
        ji2=ji2-17
        j=ipr2(ji1)+ji2
        l=2
      else
        ji1=ji1-11
        ji2=ji2-11
        l=1
        j=ipr1(ji1)+ji2
      end if
      if(k1 == 1 .and. k2 == 2) then
        if(l == 1) then
           if(i12ps(j) >= 0) then
              coef=dsqrt(dble(i12ps(j))/dble(i12pv(j)))
           else
              coef=-dsqrt(-dble(i12ps(j))/dble(i12pv(j)))
           end if
        else
           if(i12as(j) >= 0) then
              coef=dsqrt(dble(i12as(j))/dble(i12av(j)))
           else
              coef=-dsqrt(-dble(i12as(j))/dble(i12av(j)))
           end if
        end if
      else if(k1 == 0 .and. k2 == 3) then
        if(l == 1) then
           if(i03ps(j) >= 0) then
              coef=dsqrt(dble(i03ps(j))/dble(i03pv(j)))
           else
              coef=-dsqrt(-dble(i03ps(j))/dble(i03pv(j)))
           end if
        else
           if(i03as(j) >= 0) then
              coef=dsqrt(dble(i03as(j))/dble(i03av(j)))
           else
              coef=-dsqrt(-dble(i03as(j))/dble(i03av(j)))
           end if
        end if
      else if(k1 == 1 .and. k2 == 4) then
        if(l == 1) then
           if(i14ps(j) >= 0) then
              coef=dsqrt(dble(i14ps(j))/dble(i14pv(j)))
           else
              coef=-dsqrt(-dble(i14ps(j))/dble(i14pv(j)))
           end if
        else
           if(i14as(j) >= 0) then
              coef=dsqrt(dble(i14as(j))/dble(i14av(j)))
           else
              coef=-dsqrt(-dble(i14as(j))/dble(i14av(j)))
           end if
        end if
      else if(k1 == 0 .and. k2 == 5) then
        if(l == 1) then
           if(i05ps(j) >= 0) then
              coef=dsqrt(dble(i05ps(j))/dble(i05pv(j)))
           else
              coef=-dsqrt(-dble(i05ps(j))/dble(i05pv(j)))
           end if
        else
           if(i05as(j) >= 0) then
              coef=dsqrt(dble(i05as(j))/dble(i05av(j)))
           else
              coef=-dsqrt(-dble(i05as(j))/dble(i05av(j)))
           end if
        end if
      else if(k1 == 1 .and. k2 == 6) then
        if(l == 1) then
           if(i16ps(j) >= 0) then
              coef=dsqrt(dble(i16ps(j))/dble(i16pv(j)))
           else
              coef=-dsqrt(-dble(i16ps(j))/dble(i16pv(j)))
           end if
        else
           if(i16as(j) >= 0) then
              coef=dsqrt(dble(i16as(j))/dble(i16av(j)))
           else
              coef=-dsqrt(-dble(i16as(j))/dble(i16av(j)))
           end if
        end if
      else if(k1 == 0 .and. k2 == 7) then
        if(l == 1) then
           if(i07ps(j) >= 0) then
              coef=dsqrt(dble(i07ps(j))/dble(i07pv(j)))
           else
              coef=-dsqrt(-dble(i07ps(j))/dble(i07pv(j)))
           end if
        else
           if(i07as(j) >= 0) then
              coef=dsqrt(dble(i07as(j))/dble(i07av(j)))
           else
              coef=-dsqrt(-dble(i07as(j))/dble(i07av(j)))
           end if
        end if
      else
        write(0,"(A,4I5)") " J1 J2 = ",j1,j2
        write(0,"(A)") " ERROR IN SUB. RMEW7JJ "
        stop
      end if
      if(mod(ifaz,4) /= 0)coef=-coef
      return
      end subroutine rmew7bjj
end module rmew7bjj_i

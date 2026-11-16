!*******************************************************************
!                                                                  *
      integer function nmtejj(i2q,i2j,j,nk,nd)
!                                                                  *
!     ------------  SECTION METWO    SUBPROGRAM 21  -------------  *
!                                                                  *
!     NO FUNCTION CALLED                                           *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mtjj_c
      use mtjj2_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: i2q, i2j, j, nk, nd
!      DIMENSION LP(9),LG(9),LP3(27),LG3(27)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                :: i2v, inn, il, jp, jg, jj
      integer  :: lp(9), lg(9)
      integer :: lp3(27), lg3(27)
!-----------------------------------------------
      data lp/1,0,3,0,6,0,12,0,26/
      data lg/2,0,5,0,11,0,25,0,63/
      data lp3/1,0,8,0,16,0,25,0,35,0,46,0,58,0,71,0,85,0,100,0, &
      116,0,133,0,151,0,170/
      data lg3/7,0,15,0,24,0,34,0,45,0,57,0,70,0,84,0,99,0,115,0,&
      132,0,150,0,169,0,189/
      nmtejj=0
      if(j > 37)return
      i2v=j+1-2*i2q
      inn=(i2q*100+i2v)*100+i2j
      if(j < 9) then
         jp=lp(j)
         if(jp == 0)return
         jg=lg(j)
         if(jg == 0)return
         jj=jp
    2    if (inn == mt(jj)) go to 1
    3    jj=jj+1
         if(jj > jg)return
         go to 2
    1    nmtejj=jj
      elseif(j == 9) then
        if(max0(nk,nd) < 3) then
          jp=1
          jg=6
          jj=jp
    6     if (inn == mt9(jj)) go to 5
    4     jj=jj+1
          if(jj > jg)return
          go to 6
    5     nmtejj=jj+300
        else
          print*, "ERROR in FUNCTION NMTEJJ"
          stop
        end if
      else
        il=j-10
        jp=lp3(il)
        jg=lg3(il)
        jj=jp
   22   if (inn == mt11(jj)) go to 21
   23   jj=jj+1
        if(jj > jg)return
        go to 22
   21   nmtejj=jj
      endif
      return
      end function nmtejj

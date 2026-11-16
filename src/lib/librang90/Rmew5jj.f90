!*******************************************************************
!                                                                  *
      subroutine rmew5jj(j1,j2,k1,k2,coef)
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
      use ribojj_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: j1, j2, k1, k2
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                 :: ji1, ji2
      integer   :: i00(6), i01(6)
      integer :: i12p(3, 3),i12a(3, 3),i03p(3, 3), &
                       i03a(3, 3),i14p(3, 3),i14a(3, 3),i05p(3, 3),i05a(3, 3)
!-----------------------------------------------
      data i00/54,12,30,12,30,54/
      data i01/126,12,198,0,48,288/
      data i12p/420,360,270,360,2*0,-270,2*0/
      data i12a/0,1960,0,-1960,-1000,2430,0,2430,1980/
      data i03p/-882,3*0,384,-400,0,400,286/
      data i03a/4*0,162,-300,0,-300,22/
      data i14p/3780,-720,-4950,-720,2*0,4950,2*0/
      data i14a/2*0,3528,0,2430,1980,-3528,1980,-7722/
      data i05p/-1386,4*0,-440,0,440,-1430/
      data i05a/5*0,330,0,330,572/
!
      coef=zero
      if(imptjj(j1) /= imptjj(j2)) return
      if(j1 < 6 .or. j1 > 11) return
      if(k1 == 0 .and. k2 == 0) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i00(j1-5)))
      elseif(k1 == 1 .and. k2 == 0) then
        if(j1 /= j2) return
        if(j1 == 6) then
           coef=-dsqrt(dble(48))
        elseif(j1 == 9) then
           coef=-dsqrt(dble(20))
        elseif(j1 == 10) then
           coef=-dsqrt(dble(10))
        elseif(j1 == 11) then
           coef=-dsqrt(dble(18))
        endif
      elseif(k1 == 0 .and. k2 == 1) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i01(j1-5))/dble(7))
      elseif(k1 == 1 .and. k2 == 2) then
        if(j1 < 9) then
           ji1=j1-5
           ji2=j2-5
           if(i12p(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i12p(ji1,ji2))/dble(7))
           else
              coef=-dsqrt(-dble(i12p(ji1,ji2))/dble(7))
           endif
        else
           ji1=j1-8
           ji2=j2-8
           if(i12a(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i12a(ji1,ji2))/dble(49))
           else
              coef=-dsqrt(-dble(i12a(ji1,ji2))/dble(49))
           endif
        endif
      elseif(k1 == 0 .and. k2 == 3) then
        if(j1 < 9) then
           ji1=j1-5
           ji2=j2-5
           if(i03p(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i03p(ji1,ji2))/dble(21))
           else
              coef=-dsqrt(-dble(i03p(ji1,ji2))/dble(21))
           endif
        else
           ji1=j1-8
           ji2=j2-8
           if(i03a(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i03a(ji1,ji2))/dble(7))
           else
              coef=-dsqrt(-dble(i03a(ji1,ji2))/dble(7))
           endif
        endif
      elseif(k1 == 1 .and. k2 == 4) then
        if(j1 < 9) then
           ji1=j1-5
           ji2=j2-5
           if(i14p(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i14p(ji1,ji2))/dble(35))
           else
              coef=-dsqrt(-dble(i14p(ji1,ji2))/dble(35))
           endif
        else
           ji1=j1-8
           ji2=j2-8
           if(i14a(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i14a(ji1,ji2))/dble(49))
           else
              coef=-dsqrt(-dble(i14a(ji1,ji2))/dble(49))
           endif
        endif
      elseif(k1 == 0 .and. k2 == 5) then
        if(j1 < 9) then
           ji1=j1-5
           ji2=j2-5
           if(i05p(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i05p(ji1,ji2))/dble(21))
           else
              coef=-dsqrt(-dble(i05p(ji1,ji2))/dble(21))
           endif
        else
           ji1=j1-8
           ji2=j2-8
           if(i05a(ji1,ji2) >= 0) then
              coef=dsqrt(dble(i05a(ji1,ji2))/dble(7))
           else
              coef=-dsqrt(-dble(i05a(ji1,ji2))/dble(7))
           endif
        endif
      else
         write(0,'(A,4I5)') ' J1 J2 = ',j1,j2
         write(0,'(A)') ' ERROR IN SUB. RMEW5JJ '
         stop
      endif
      return
      end subroutine rmew5jj

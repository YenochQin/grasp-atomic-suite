!*******************************************************************
!                                                                  *
      subroutine rmew3jj(j1,j2,k1,k2,coef)
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
      integer               :: ji1
      integer :: i00(3), i10(3), i01(3)
!-----------------------------------------------
      data i00/16,6,10/
      data i10/12,12,0/
      data i01/12,0,12/
!
      coef=zero
      if(imptjj(j1) /= imptjj(j2)) return
      if(j1 < 2 .or. j1 > 5) return
      ji1=j1-2
      if(k1 == 0 .and. k2 == 0) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i00(ji1)))
      elseif(k1 == 1 .and. k2 == 0) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i10(ji1)))
      elseif(k1 == 0 .and. k2 == 1) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i01(ji1)))
      elseif(k1 == 1 .and. k2 == 2) then
        if(j1 == 3 .and. j2 == 3) then
           coef=dsqrt(dble(60))
        elseif(j1 == 5 .and. j2 == 4) then
           coef=dsqrt(dble(30))
        elseif(j1 == 4 .and. j2 == 5) then
           coef=-dsqrt(dble(30))
        endif
      elseif(k1 == 0 .and. k2 == 3) then
        if(j1 == 3 .and. j2 == 3) then
           coef=-dsqrt(dble(28))
        elseif(j1 == 5 .and. j2 == 5) then
           coef=dsqrt(dble(28))
        endif
      else
        write(0,'(A,4I5)') ' J1 J2 = ',j1,j2
        write(0,'(A)') ' ERROR IN SUB. RMEW3JJ '
        stop
      endif
      return
      end subroutine rmew3jj

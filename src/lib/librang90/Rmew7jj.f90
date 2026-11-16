module rmew7jj_i
    implicit none (type, external)
    private
    public :: rmew7jj
contains
!*******************************************************************
!                                                                  *
      subroutine rmew7jj(j1,j2,k1,k2,coef)
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
      use rmew7bjj_i, only: rmew7bjj
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: j1, j2, k1, k2
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i00(14), i10(14), i01(14)
!-----------------------------------------------
      data i00/32,48,128,80,96,128,20,60,20,108,36,44,156,68/
      data i10/6,9,120,15,18,24,2*30,0,54,2*0,78,0/
      data i01/10,35,168,165,286,680,0,30,10,180,60,110,546,408/
      coef=zero
      if(imptjj(j1) /= imptjj(j2)) return
      if(j1 < 12 .or. j1 > 25) return
      if(k1 == 0 .and. k2 == 0) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i00(j1-11)))
      else if(k1 == 1 .and. k2 == 0) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i10(j1-11)))
      else if(k1 == 0 .and. k2 == 1) then
        if(j1 /= j2) return
        coef=-dsqrt(dble(i01(j1-11))/dble(7))
      else
        call rmew7bjj(j1,j2,k1,k2,coef)
      end if
      return
      end subroutine rmew7jj
end module rmew7jj_i

module rmew1jj_i
    implicit none (type, external)
    private
    public :: rmew1jj
contains
!*******************************************************************
!                                                                  *
      subroutine rmew1jj(j1,j2,k1,k2,coef)
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
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: j1, j2, k1, k2
      real(real64), intent(out) :: coef
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i10(2), i01(2)
!-----------------------------------------------
      data i01/6,0/
      data i10/0,6/
!
      coef=zero
      if(imptjj(j1) /= imptjj(j2)) return
      if(j1 > 2) return
      if(k1 == 0 .and. k2 == 0) then
        coef=-dsqrt(dble(2))
      else if(k1 == 1 .and. k2 == 0) then
        coef=-dsqrt(dble(i10(j1)))
      else if(k1 == 0 .and. k2 == 1) then
        coef=-dsqrt(dble(i01(j1)))
      else
        write(0,"(A,4I5)") " J1 J2 = ",j1,j2
        write(0,"(A)") " ERROR IN SUB. RMEW1JJ "
        stop
      end if
      return
      end subroutine rmew1jj
end module rmew1jj_i

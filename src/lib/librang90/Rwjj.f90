!*******************************************************************
!                                                                  *
      subroutine rwjj(j,j1,j2,k1,k2,coef)
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
      use cons_c
      use ribojj_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rmew1jj_i
      use rmew3jj_i
      use rmew5jj_i
      use rmew7jj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)  :: j, j1, j2, k1, k2
      real(real64), intent(out) :: coef
!-----------------------------------------------
      if(j == 1) then
         call rmew1jj(j1,j2,k1,k2,coef)
      elseif(j == 3) then
         call rmew3jj(j1,j2,k1,k2,coef)
      elseif(j == 5) then
         call rmew5jj(j1,j2,k1,k2,coef)
      elseif(j == 7) then
         call rmew7jj(j1,j2,k1,k2,coef)
!GG        ELSEIF(J.EQ.9) THEN
!GG          CALL SUWJJ(K1,K2,J,J1,J2,COEF)
      else
         write(0,'(A,I5)') ' KLAIDA SUB. RWJJ J=',j
         stop
      endif
      return
      end subroutine rwjj

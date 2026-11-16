!***********************************************************************
!                                                                      *
      real(real64) function fzalf (n, kappa, z)
!                                                                      *
!   An estimate of the function  F (Z*\alpha) is computed here.        *
!                                                                      *
!   Call(s) to: [RCI92]: KLAMAQ, MOHR.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1990   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mohr_i, only: mohr
      use klamaq_i, only: klamaq
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: kappa
      real(real64)  :: z
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: neff
      real(real64) :: value
!-----------------------------------------------
!
      if (n <= 2) then
         call mohr (n, kappa, z, value)
      else
         if (kappa==(-1) .or. kappa==1 .or. kappa==(-2)) then
            neff = 2
            call mohr (neff, kappa, z, value)
         else
            call klamaq (n, kappa, z, value)
         end if
      end if
!
      fzalf = value
!
      return
      end function fzalf

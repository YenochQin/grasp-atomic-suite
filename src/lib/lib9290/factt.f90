module factt_i
    implicit none (type, external)
    private
    public :: factt
contains
!***********************************************************************
!                                                                      *
      subroutine factt
!                                                                      *
!   Calculates the logs  of factorials required by the Racah coeffi-   *
!   cient routine DRACAH.                                              *
!                                                                      *
!   Written by N S Scott                    Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use facts_c, only: double, gam, mfact
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: x
!-----------------------------------------------
!
      gam(1) = 1.0_real64
      gam(2) = 1.0_real64
      x = 2.0_real64
!
      do i = 3, 30
         gam(i) = gam(i-1)*x
         x = x + 1.0_real64
      end do
!
      do i = 1, 30
         gam(i) = log(gam(i))
      end do
!
      x = 3.0e1_real64
!
      do i = 31, mfact
         gam(i) = gam(i-1) + log(x)
         x = x + 1.0_real64
      end do
!
      return
      end subroutine factt
end module factt_i

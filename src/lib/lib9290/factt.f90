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
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:26   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use facts_c
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: x
!-----------------------------------------------
!
      gam(1) = 1.0d00
      gam(2) = 1.0d00
      x = 2.0d00
!
      do i = 3, 30
         gam(i) = gam(i-1)*x
         x = x + 1.0d00
      end do
!
      do i = 1, 30
         gam(i) = log(gam(i))
      end do
!
      x = 3.0d01
!
      do i = 31, mfact
         gam(i) = gam(i-1) + log(x)
         x = x + 1.0d00
      end do
!
      return
      end subroutine factt

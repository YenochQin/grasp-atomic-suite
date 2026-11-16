!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function rintdens (i, j)
!                                                                      *
!   The value of RINTDENS is an approximation to:                      *
!                                                                      *
!                                                                      *
!      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
!                     I     J       I     J                            *
!                                                                      *
!   Call(s) to: [SMS92]:  POLINT                                       *
!                                                                      *
!   Written by Per Jonsson             Last revision: 24 Dec 1992      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:07:11   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  11/02/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use def_c,           only: vacuum_speed_of_light_au, pi
      use grid_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use polint_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: i, j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, ii, nmin
      real(double), dimension(3) :: xa, ya
      real(double) :: dens
!-----------------------------------------------
!
      nmin = 10

      ii = 1
      do l = nmin,nmin+4,2
         xa(ii) = radial_grid(l)
         ya(ii) = (radial_large_components(l,i)*radial_large_components(l,j) + radial_small_components(l,i)* &
            & radial_small_components(l,j))/(4.0d00*pi*radial_grid(l)*radial_grid(l))
         ii = ii + 1
      end do
      call polint (xa,ya,dens)
      rintdens = dens
!
      return
      end function rintdens

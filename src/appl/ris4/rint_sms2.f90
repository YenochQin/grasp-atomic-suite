!***********************************************************************
!                                                                      *
      real(real64) function rint_sms2 (i, j)
!                                                                      *
!   The value of RINT_SMS2 is an approximation to:                     *
!                                                                      *
!                                                                      *
!   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
!   Range.                                                             *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92]: QUAD.                                         *
!                                                                      *
!   Written by  G. Gaigalas                                            *
!           and E. Gaudamauskas                                        *
!                        Last             revision:  09 October 2009   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use tatb_c
      use wave_c
      use def_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use sigma_1_i
      use sigma_2_i
      use quad_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer ,intent(in)  :: i, j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: l
      real(real64) :: result, apart1, apart2
!-----------------------------------------------
!
!   Tabulate integrand as required for SUBROUTINE QUAD
!
      mtp = min (orbital_last_grid_index(i),orbital_last_grid_index(j))
!
!   Value at first tabulation point is arbitrary
!
      call sigma_1(2,i,j,apart1)
      call sigma_2(2,i,j,apart2)
      ta(1) = 0.0_real64
      do 1 l = 2,mtp
!HFS        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
        ta(l) =  radial_grid_derivative(l)*(apart1*radial_small_components(l,i)*radial_large_components(l,j) &
           & &
              -         apart2*radial_large_components(l,i)*radial_small_components(l,j))/radial_grid(l)
    1 continue
!
!   Perform integration
!
      call quad (result)
      rint_sms2 = -result*nuclear_charge/speed_of_light_au
!
      return
      end function rint_sms2

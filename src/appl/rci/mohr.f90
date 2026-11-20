!***********************************************************************
!                                                                      *
      subroutine mohr(n, kappa, z, fzalfa)
!                                                                      *
!   The  function  F (Z*alpha)  for the  1s  2s  2p-  2p  symmetries   *
!   is computed here.    A value is obtained by interpolating in, or   *
!   extrapolating from, the table due to  P J Mohr.   See  P J Mohr,   *
!   At Data Nucl Data Tables 29 (1983) 453.                            *
!                                                                      *
!   Call(s) to: [LIB92]: INTERP.                                       *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use aitken_interpolation, only: interp
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: kappa
      real(real64)  :: z
      real(real64) , intent(out) :: fzalfa
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!   Number of data points
      integer, parameter :: numval = 12
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: val1s(numval), val2s(numval), val2p1(numval), val2p3(numval), arg(numval)
      real(real64) :: accy, value
!-----------------------------------------------
!
!
!   1s data :
!
      data val1s/ 10.3168_real64, 4.6540_real64, 3.2460_real64, 2.5519_real64, 2.1351_real64, &
         1.8644_real64, 1.6838_real64, 1.5675_real64, 1.5032_real64, 1.4880_real64, 1.5317_real64, &
         1.6614_real64/
!
!   2s data:
!
      data val2s/ 10.5468_real64, 4.8930_real64, 3.5063_real64, 2.8391_real64, 2.4550_real64, &
         2.2244_real64, 2.0948_real64, 2.0435_real64, 2.0650_real64, 2.1690_real64, 2.3870_real64, &
         2.7980_real64/
!
!   2p- data:
!
      data val2p1/ -0.1264_real64, -0.1145_real64, -0.0922_real64, -0.0641_real64, -0.0308_real64, &
         0.0082_real64, 0.0549_real64, 0.1129_real64, 0.1884_real64, 0.2934_real64, 0.4530_real64, &
         0.7250_real64/
!
!   2p data:
!
      data val2p3/ 0.1235_real64, 0.1303_real64, 0.1436_real64, 0.1604_real64, 0.1794_real64, &
         0.1999_real64, 0.2215_real64, 0.2440_real64, 0.2671_real64, 0.2906_real64, 0.3141_real64, &
         0.3367_real64/
!
!   Z data:
!
      data arg/ 1.0_real64, 10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64, 60.0_real64, &
         70.0_real64, 80.0_real64, 90.0_real64, 100.0_real64, 110.0_real64/
!
!----------------------------------------------------------------------*
!
!   Convergence criterion for interpolation
!
      data accy/ 1.0e-3_real64/
!
!   Interpolate or issue error message as appropriate
!
      if (n == 1) then
         if (kappa == (-1)) then
            call interp (arg, val1s, numval, z, value, accy)
         else
            write (*, '("MOHR:")')
            write (*, '(" Principal quantum number, ",i12,", kappa, ",1i3,".")') n, kappa
            stop
         end if
      else if (n == 2) then
         select case (kappa)
         case (-1)
            call interp (arg, val2s, numval, z, value, accy)
         case (1)
            call interp (arg, val2p1, numval, z, value, accy)
         case (-2)
            call interp (arg, val2p3, numval, z, value, accy)
         case default
            write (*, '("MOHR:")')
            write (*, '(" Principal quantum number, ",i12,", kappa, ",1i3,".")') n, kappa
            stop
         end select
      else
         write (*, '("MOHR:")')
         write (*, '(" Principal quantum number, ",1i2,", Should be either 1 or 2.")') n
         stop
      end if
!
      fzalfa = value
!
      return
!
      return
!
      end subroutine mohr

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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use interp_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: kappa
      real(double)  :: z
      real(double) , intent(out) :: fzalfa
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!   Number of data points
      integer, parameter :: numval = 12
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(numval) :: val1s, val2s, val2p1, val2p3, arg
      real(double) :: accy, value
!-----------------------------------------------
!
!
!   1s data :
!
      data val1s/ 10.3168d00, 4.6540d00, 3.2460d00, 2.5519d00, 2.1351d00, &
         1.8644d00, 1.6838d00, 1.5675d00, 1.5032d00, 1.4880d00, 1.5317d00, &
         1.6614d00/
!
!   2s data:
!
      data val2s/ 10.5468d00, 4.8930d00, 3.5063d00, 2.8391d00, 2.4550d00, &
         2.2244d00, 2.0948d00, 2.0435d00, 2.0650d00, 2.1690d00, 2.3870d00, &
         2.7980d00/
!
!   2p- data:
!
      data val2p1/ -0.1264d00, -0.1145d00, -0.0922d00, -0.0641d00, -0.0308d00, &
         0.0082d00, 0.0549d00, 0.1129d00, 0.1884d00, 0.2934d00, 0.4530d00, &
         0.7250d00/
!
!   2p data:
!
      data val2p3/ 0.1235d00, 0.1303d00, 0.1436d00, 0.1604d00, 0.1794d00, &
         0.1999d00, 0.2215d00, 0.2440d00, 0.2671d00, 0.2906d00, 0.3141d00, &
         0.3367d00/
!
!   Z data:
!
      data arg/ 1.0d00, 10.0d00, 20.0d00, 30.0d00, 40.0d00, 50.0d00, 60.0d00, &
         70.0d00, 80.0d00, 90.0d00, 100.0d00, 110.0d00/
!
!----------------------------------------------------------------------*
!
!   Convergence criterion for interpolation
!
      data accy/ 1.0d-03/
!
!   Interpolate or issue error message as appropriate
!
      if (n == 1) then
         if (kappa == (-1)) then
            call interp (arg, val1s, numval, z, value, accy)
         else
            write (*, 300)
            write (*, 301) n, kappa
            stop
         endif
      else if (n == 2) then
         select case (kappa)
         case (-1)
            call interp (arg, val2s, numval, z, value, accy)
         case (1)
            call interp (arg, val2p1, numval, z, value, accy)
         case (-2)
            call interp (arg, val2p3, numval, z, value, accy)
         case default
            write (*, 300)
            write (*, 301) n, kappa
            stop
         end select
      else
         write (*, 300)
         write (*, 302) n
         stop
      endif
!
      fzalfa = value
!
      return
!
  300 format('MOHR:')
  301 format(' Principal quantum number, ',i12,', kappa, ',1i3,'.')
  302 format(' Principal quantum number, ',1i2,', Should be either 1 or 2.')
      return
!
      end subroutine mohr

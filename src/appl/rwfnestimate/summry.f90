

!***********************************************************************
!                                                                      *
      subroutine summry(nunit)
!                                                                      *
!   Prints  a summary of the complete  list of subshell radial wave-   *
!   functions on NUNIT.                                                *
!                                                                      *
!   Call(s) to: [LIB92]:                                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!...Modified by Yanting Li, output <r> for each orbital, 19 Aug 2022
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use orb_c
      use wave_c
      use whfrom_c, only: source
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rint_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nunit
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, lenth
!-----------------------------------------------
!
!
      write (nunit, 300)
!
      do i = 1, nw
         lenth = len_trim(source(i))
!--------Output P(2) and Q(2) for each orbital
!         WRITE (NUNIT, 301) NP(I), NH(I), E(I), PZ(I), GAMA(I), PF(2,I), QF(2,I&
!            ), MF(I), SOURCE(I)(1:LENTH)
!--------Output <r> (RINT(I,I,1)) for each orbital
         write (nunit, 301) np(i), nh(i), e(i), pz(i), gama(i), rint(i,i,1)&
             , mf(i), source(i)(1:lenth)
!         WRITE (NUNIT,302) SOURCE(I)(1:LENTH)
      end do
!
      return
!
!-Format for P(2) and Q(2)
!  300 FORMAT('Shell',6X,'e',11X,'p0',8X,'gamma',8X,'P(2)',7X,'Q(2)',6X,'MTP',&
!         '  SRC'/)
!  301 FORMAT(1X,I2,A2,5D12.4,I5,2X,A3)
!-Format for <r>
  300 format('Shell',6x,'e',11x,'p0',8x,'gamma',8x,'<r>',6x,'MTP','  SRC'/)
  301 format(1x,i2,a2,4d12.4,i5,2x,a3)
      return
!  302 FORMAT ('      Source: ',A)
!
      end subroutine summry

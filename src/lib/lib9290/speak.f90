!***********************************************************************
!                                                                      *
      subroutine speak(ja, jb, ia1, ib1, ia2, ib2, k, x)
!                                                                      *
!   Output MCP  coefficients and integral parameters to COMMON block   *
!   /BUFFER/. Also print these if  IBUG1 = 1 .                         *
!                                                                      *
!                                           Last Update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:59:36   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use buffer_c,        only: label, coeff, nbdim, nvcoef
      use debug_c ,        only: ibug1
      use orb_c,           only: np, nh
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja
      integer, intent(in) :: jb
      integer, intent(in) :: ia1
      integer, intent(in) :: ib1
      integer, intent(in) :: ia2
      integer, intent(in) :: ib2
      integer, intent(in) :: k
      real(double), intent(in) :: x
!-----------------------------------------------
!
!
      if (ibug1 /= 0) write (99, 300) ja, jb, np(ia1), nh(ia1), np(ib1), nh(ib1&
         ), np(ia2), nh(ia2), np(ib2), nh(ib2), k, x
!
!   Increment counter
!
      nvcoef = nvcoef + 1
!
!   Ensure that arrays are of adequate size; reallocate if necessary
!
      if (nvcoef > nbdim) call alcbuf (2)
!
!   Store integral indices and coefficient in COMMON/BUFFER/
!
      label(1,nvcoef) = ia1
      label(2,nvcoef) = ib1
      label(3,nvcoef) = ia2
      label(4,nvcoef) = ib2
      label(5,nvcoef) = k
      coeff(nvcoef) = x
!
      return
!
  300 format(2(1x,1i2),4(1x,i2,a2),1x,i2,1x,1p,d19.12)
      return
!
      end subroutine speak

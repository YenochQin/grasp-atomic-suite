!***********************************************************************
!                                                                      *
      subroutine talk(ja, jb, nu, ia, ib, ic, id, itype, coef)
!                                                                      *
!   Print  coefficients  and  integral  parameters  if IBUG1 > 0 and   *
!   write to disk.                                                     *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:16:29   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only:  keyorb
      use buffer_c
      use debug_c
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use alcbuf_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ja, jb, nu, ia, ib, ic, id, itype
      real(double), intent(in) :: coef
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = keyorb
!-----------------------------------------------
!   Print coefficient if requested
!
      if (ibug1 /= 0) write (99, 300) ja, jb, np(ia), nh(ia), np(ib), nh(ib), &
         np(ic), nh(ic), np(id), nh(id), nu, itype, coef
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
      label(1,nvcoef) = ia
      label(2,nvcoef) = ib
      label(3,nvcoef) = ic
      label(4,nvcoef) = id
      label(5,nvcoef) = nu
      label(6,nvcoef) = itype
      coeff(nvcoef) = coef
!
      return
!
  300 format(2(1x,1i2),4(1x,i2,a2),1x,1i2,1x,1i2,1x,1p,d19.12)
      return
!
      end subroutine talk

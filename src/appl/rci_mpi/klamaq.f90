!***********************************************************************
!                                                                      *
      subroutine klamaq(n, kappa, z, fzalfa)
!                                                                      *
!   The function  F (Z*\alpha)  is estimated here. We use the series   *
!   expansion given by  Eqs (1) and (2) and the table of Bethe loga-   *
!   rithms in  S Klarsfeld and A Maquet, Physics Letters  43B (1973)   *
!   201. The vacuum-polarization contribution in Eq (2) is omitted.    *
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
      use vast_kind_param, only: double
      use def_c,           only: speed_of_light_au
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: kappa
      real(double) , intent(in) :: z
      real(double) , intent(out) :: fzalfa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, loc
      real(double), dimension(36) :: bethe
      real(double) :: c401, c402, ovlfac, bethel, term, zalfa, factor
      logical :: first
!-----------------------------------------------
!
      data bethe/ 2.9841285d00, 2.8117699d00, -0.0300167d00, 2.7676636d00, &
         -0.0381902d00, -0.0052321d00,  2.7498118d00, -0.0419549d00, &
         -0.0067409d00, -0.0017337d00,  2.7408237d00, -0.0440347d00, &
         -0.0076008d00, -0.0022022d00, -0.0007721d00,  2.7356642d00, &
         -0.0453122d00, -0.0081472d00, -0.0025022d00, -0.0009628d00, &
         -0.0004079d00,  2.7324291d00, -0.0461552d00, -0.0085192d00, &
         -0.0027091d00, -0.0010945d00, -0.0004997d00, -0.0002409d00, &
          2.7302673d00, -0.0467413d00, -0.0087850d00, -0.0028591d00, &
         -0.0011904d00, -0.0005665d00, -0.0002904d00, -0.0001539d00/
!
!----------------------------------------------------------------------*
!
      data first/ .true./
!
      data c401/ 0.0d00/
      data c402/ 0.0d00/
      data ovlfac/ 0.0d00/
!
!   Set up the constants
!
      if (first) then
!
         c401 = 11.0d00/24.0d00
         c402 = 3.0d00/8.0d00
         ovlfac = 4.0d00/3.0d00
!
         first = .false.
!
      endif
!
!   Ensure that the principal quantum number is in range
!
      if (n<1 .or. n>8) then
         write (*, 300)
         write (*, 301) n
         stop
      endif
!
!   Determine the azimuthal quantum number
!
      if (kappa > 0) then
         l = kappa
      else if (kappa == 0) then
         write (*, 300)
         write (*, 302)
         stop
      else
         l = (-kappa) - 1
      endif
!
!   Ensure that the azimuthal quantum number is in range
!
      if (l > n - 1) then
         write (*, 300)
         write (*, 303) kappa, n
         stop
      endif
!
!   Find the appropriate entry in the table
!
      loc = (n*n - n)/2 + l + 1
      bethel = bethe(loc)
!
!   Determine the quantity in square brackets in eq. (1) of
!   Klarsfeld and Maquet
!
      term = -bethel
!
      if (kappa > 0) then
         term = term - c402/dble(l*(l + l + 1))
      else
         term = term + c402/dble((l + 1)*(l + l + 1))
         if (kappa == (-1)) then
            zalfa = z/speed_of_light_au
            factor = dlog(zalfa)
            factor = -(factor + factor)
            term = term + factor + c401
         endif
      endif
!
      fzalfa = ovlfac*term
!
      return
!
  300 format('KLAMAQ:')
  301 format(' Principal quantum number, ',1i2,&
         ', should be in the range  1--8.')
  302 format(' Kappa is  0 .')
  303 format(' Kappa, ',1i3,', is out of range for n, ',1i2,'.')
      return
!
      end subroutine klamaq

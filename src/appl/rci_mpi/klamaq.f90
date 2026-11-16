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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use def_c,           only: speed_of_light_au
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n
      integer , intent(in) :: kappa
      real(real64) , intent(in) :: z
      real(real64) , intent(out) :: fzalfa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, loc
      real(real64) :: bethe(36)
      real(real64) :: c401, c402, ovlfac, bethel, term, zalfa, factor
      logical :: first
!-----------------------------------------------
!
      data bethe/ 2.9841285_real64, 2.8117699_real64, -0.0300167_real64, 2.7676636_real64, &
         -0.0381902_real64, -0.0052321_real64,  2.7498118_real64, -0.0419549_real64, &
         -0.0067409_real64, -0.0017337_real64,  2.7408237_real64, -0.0440347_real64, &
         -0.0076008_real64, -0.0022022_real64, -0.0007721_real64,  2.7356642_real64, &
         -0.0453122_real64, -0.0081472_real64, -0.0025022_real64, -0.0009628_real64, &
         -0.0004079_real64,  2.7324291_real64, -0.0461552_real64, -0.0085192_real64, &
         -0.0027091_real64, -0.0010945_real64, -0.0004997_real64, -0.0002409_real64, &
          2.7302673_real64, -0.0467413_real64, -0.0087850_real64, -0.0028591_real64, &
         -0.0011904_real64, -0.0005665_real64, -0.0002904_real64, -0.0001539_real64/
!
!----------------------------------------------------------------------*
!
      data first/ .true./
!
      data c401/ 0.0_real64/
      data c402/ 0.0_real64/
      data ovlfac/ 0.0_real64/
!
!   Set up the constants
!
      if (first) then
!
         c401 = 11.0_real64/24.0_real64
         c402 = 3.0_real64/8.0_real64
         ovlfac = 4.0_real64/3.0_real64
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

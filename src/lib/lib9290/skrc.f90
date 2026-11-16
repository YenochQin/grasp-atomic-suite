module skrc_i
    implicit none (type, external)
    private
    public :: skrc
contains
!***********************************************************************
!                                                                      *
      subroutine skrc(is, kaps, ks, kd1, kd2, ke1, ke2)
!                                                                      *
!   Determines the range of the tensor rank k for Coulomb integral.    *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: kd1
      integer, intent(out) :: kd2
      integer, intent(out) :: ke1
      integer, intent(out) :: ke2
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      integer, intent(in) :: ks(4)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isd1, isd2, kd1a, kd1b, kd2a, kd2b, ise1, ise2, ke1a, &
         ke1b, ke2a, ke2b
!-----------------------------------------------
!
!
      kd2 = 0
      ke2 = 0
!
!   Direct terms --- KD1 = minimum k , KD2 = number of terms
!
      isd1 = 1
      if (kaps(1)*kaps(3) < 0) isd1 = -1
      isd2 = 1
      if (kaps(2)*kaps(4) < 0) isd2 = -1
      kd1a = abs(ks(1)-ks(3))
      if (isd1 < 0) kd1a = kd1a + 2
      kd1b = abs(ks(2)-ks(4))
      if (isd2 < 0) kd1b = kd1b + 2
      if (mod((kd1a - kd1b)/2,2) == 0) then
         kd2a = ks(1) + ks(3) - 2
         if (isd1 > 0) kd2a = kd2a - 2
         kd2b = ks(2) + ks(4) - 2
         if (isd2 > 0) kd2b = kd2b - 2
         kd1 = max(kd1a,kd1b)/2
         kd2 = min(kd2a,kd2b)/2
         kd2 = (kd2 - kd1)/2 + 1
      end if
!
!   Exchange terms --- KE1 = minimum k , KE2 = number of terms
!
      if (is(1)==is(2) .or. is(3)==is(4)) return
      ise1 = 1
      if (kaps(1)*kaps(4) < 0) ise1 = -1
      ise2 = 1
      if (kaps(2)*kaps(3) < 0) ise2 = -1
      ke1a = abs(ks(1)-ks(4))
      if (ise1 < 0) ke1a = ke1a + 2
      ke1b = abs(ks(2)-ks(3))
      if (ise2 < 0) ke1b = ke1b + 2
      if (mod((ke1a - ke1b)/2,2) /= 0) return
      ke2a = ks(1) + ks(4) - 2
      if (ise1 > 0) ke2a = ke2a - 2
      ke2b = ks(2) + ks(3) - 2
      if (ise2 > 0) ke2b = ke2b - 2
      ke1 = max(ke1a,ke1b)/2
      ke2 = min(ke2a,ke2b)/2
      ke2 = (ke2 - ke1)/2 + 1
!
      return
      end subroutine skrc
end module skrc_i

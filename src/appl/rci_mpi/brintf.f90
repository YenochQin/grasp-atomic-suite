!***********************************************************************
!                                                                      *
      real(real64) function brintf (itype, ia, ib, ic, id, k)
!                                                                      *
!   Computes integrals for the transverse photon interaction.          *
!                                                                      *
!   Call(s) to: [RCI92]: BESSEL, BRRA.                                 *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def,   only: orbital_index_radix
      use stor_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use bessel_i
      use brra_i
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: itype
      integer, intent(in) :: ia
      integer, intent(in) :: ib
      integer, intent(in) :: ic
      integer, intent(in) :: id
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: key = orbital_index_radix
      integer, parameter :: key2 = key*key
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icod, icod1, icod2, icod3, icod4
!-----------------------------------------------
!
      go to (1,4,3,6,2,5) itype
!
!   Type 1 and 5 integrals require j(k), n(k) Bessel fuctions
!   Type 5 integrals only require w = wab Bessel functions
!
    1 continue
      if (ia/=ic .or. ib/=id .or. ia/=id .or. ic==ib) call bessel (ic, id, 1, 2&
         , k)
    2 continue
      call bessel (ia, ib, 1, 1, k)
      go to 6
!
!   Type 3 integrals require j(k), n(k) Bessel functions for either
!   w = wab or w = cd whichever is non-zero.
!
    3 continue
      if (ia /= ib) call bessel (ia, ib, 1, 1, k)
      if (ic /= id) call bessel (ic, id, 1, 2, k)
      go to 6
!
!   Type 2 and 6 integrals require j(k), n(k) and j(k+2), n(k+2)
!   Bessel fuctions
!   Type 6 integrals only require w = wab Bessel functions.
!
    4 continue
      if (ia/=ic .or. ib/=id .or. ia/=id .or. ic/=ib) then
!
         icod = max(ic,id) + key*min(ic,id)
         icod1 = icod + key2*(k - 1)
         icod2 = icod + key2*(k + 1)
         icod = max(ia,ib) + key*min(ia,ib)
         icod3 = icod + key2*(k - 1)
         icod4 = icod + key2*(k + 1)
         if (icod1==keep(1,2) .and. icod2==keep(2,2) .and. icod3==keep(1,1)&
             .and. icod4==keep(2,1)) go to 6
         if (icod1==keep(1,1) .and. icod2==keep(2,1) .and. icod3==keep(1,2)&
             .and. icod4==keep(2,2)) go to 6
         call bessel (ic, id, 1, 2, k - 1)
         call bessel (ic, id, 2, 2, k + 1)
      endif
!
    5 continue
      call bessel (ia, ib, 1, 1, k - 1)
      call bessel (ia, ib, 2, 1, k + 1)
!
!   Compute the integral
!
    6 continue
      brintf = brra(itype,ia,ib,ic,id,k)
!
      return
      end function brintf

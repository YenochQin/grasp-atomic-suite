module snrc_i
    implicit none (type, external)
    private
    public :: snrc
contains
!***********************************************************************
!                                                                      *
      subroutine snrc(is, kaps, ks, nd1, nd2, ne1, ne2, ibrd, ibre)
!                                                                      *
!   Determines the range of tensor rank NU for direct/exchange terms,  *
!   and classifies the types of radial integral.                       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS      : Orbital labels                                        *
!      KAPS    : Values of 2*kappa                                     *
!      KS      : Values of 2*J+1                                       *
!                                                                      *
!   Outputs:                                                           *
!                                                                      *
!      ND1/NE1 : Lowest NU value for direct/exchange types             *
!      ND2/NE2 : Corresponding number of contributing NU values: NU    *
!                = ND1, ND1+2 ,..., ND1+2*(ND2-1) etc                  *
!      IBRD    : Classify types of  radial  integrals  contributing;   *
!      IBRE      negative value implies null contribution              *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: nd1
      integer, intent(out) :: nd2
      integer, intent(out) :: ne1
      integer, intent(out) :: ne2
      integer, intent(out) :: ibrd
      integer, intent(out) :: ibre
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      integer, intent(in) :: ks(4)
!----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: iac, iad, nd1a, nd2a, ne1a, ne2a
!-----------------------------------------------
!
      nd2 = 0
      ne2 = 0
!
!   2.0  Form limits for direct terms
!
      iac = 1
      if (kaps(1)*kaps(3) < 0) iac = -1
      iad = 1
      if (kaps(2)*kaps(4) < 0) iad = -1
      nd1 = abs(ks(1)-ks(3))/2 - 1
      if (iac == (-1)) nd1 = nd1 + 1
      if (nd1 == (-1)) nd1 = 1
      nd1a = abs(ks(2)-ks(4))/2 - 1
      if (iad == (-1)) nd1a = nd1a + 1
      if (nd1a == (-1)) nd1a = 1
      if (mod(nd1 - nd1a,2) /= 0) then
         ibrd = -1
      else
         nd2 = abs(ks(1)+ks(3))/2
         if (iac == (-1)) nd2 = nd2 + 1
         nd2a = abs(ks(2)+ks(4))/2
         if (iad == (-1)) nd2a = nd2a + 1
         nd1 = max(nd1,nd1a)
         nd2 = min(nd2,nd2a)
         nd2 = (nd2 - nd1)/2 + 1
!
!   2.1  Identify type of radial integrals
!
         ibrd = 1
         if (is(1)==is(3) .and. is(2)/=is(4) .or. is(1)/=is(3) .and. is(2)==is(&
            4)) ibrd = 2
         if (is(1)==is(3) .and. is(2)==is(4)) ibrd = 3
      end if
!
!   3.0  Form limits for exchange terms
!
      if (is(1)==is(2) .or. is(3)==is(4)) then
         ibre = -1
         return
      end if
      iac = 1
      if (kaps(1)*kaps(4) < 0) iac = -1
      iad = 1
      if (kaps(2)*kaps(3) < 0) iad = -1
      ne1 = iabs(ks(1)-ks(4))/2 - 1
      if (iac == (-1)) ne1 = ne1 + 1
      if (ne1 == (-1)) ne1 = 1
      ne1a = abs(ks(2)-ks(3))/2 - 1
      if (iad == (-1)) ne1a = ne1a + 1
      if (ne1a == (-1)) ne1a = 1
      if (mod(ne1 - ne1a,2) /= 0) then
         ibre = -1
         return
      end if
!
      ne2 = abs(ks(1)+ks(4))/2
      if (iac == (-1)) ne2 = ne2 + 1
      ne2a = abs(ks(2)+ks(3))/2
      if (iad == (-1)) ne2a = ne2a + 1
      ne1 = max(ne1,ne1a)
      ne2 = min(ne2,ne2a)
      ne2 = (ne2 - ne1)/2 + 1
!
!   3.1  Identify type of radial integrals
!
      ibre = 1
      if (is(1)==is(4) .and. is(2)/=is(3) .or. is(1)/=is(4) .and. is(2)==is(3)&
         ) ibre = 2
      if (is(1)==is(3) .and. is(2)==is(4)) ibre = 4
      return
!
      end subroutine snrc
end module snrc_i

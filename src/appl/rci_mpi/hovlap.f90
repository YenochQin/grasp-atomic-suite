!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function hovlap (p, q, mtpo, np, kappa, z)
!                                                                      *
!   This subprogram computes the overlap of the orbital tabulated in   *
!   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
!   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
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
      use parameter_def,   only: nnnp
      use grid_c
      use horb_c, only: ph, qh
      use tatb_c, only: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dcbsrw_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mtpo
      integer  :: np
      integer  :: kappa
      real(double)  :: z
      real(double), dimension(nnnp), intent(in) :: p
      real(double), dimension(nnnp), intent(in) :: q
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtph, i
      real(double) :: eh, pzh, result
!-----------------------------------------------
!
!   Set up the hydrogenic orbital
!
      call dcbsrw (np, kappa, z, eh, pzh, ph, qh, mtph)
!
!   Compute the overlap
!
      mtp = min(mtph,mtpo)
      ta(1) = 0.0d00
      ta(2:mtp) = (p(2:mtp)*ph(2:mtp)+q(2:mtp)*qh(2:mtp))*rp(2:mtp)
      call quad (result)
!
      hovlap = result
!
      return
      end function hovlap

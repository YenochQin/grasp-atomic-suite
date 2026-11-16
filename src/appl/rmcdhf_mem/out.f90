!***********************************************************************
!                                                                      *
      subroutine out(j, jp, p, q)
!                                                                      *
!   This subroutine carries out the step-by-step outward integration   *
!   of a pair of inhomogeneous Dirac radial equations.                 *
!                                                                      *
!   arguments:                                                         *
!                                                                      *
!      J:   (Input) Orbital index of function to be computed           *
!      JP:  (Input) The join point; the outward integration stops      *
!           at this tabulation index                                   *
!      P,Q: (Input and output) on input, elements 1 to 3 of both       *
!           arrays must be tabulated; on output, the arrays are        *
!           tabulated up to point JP                                   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 08 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      use grid_c,           only: transformed_grid_step, radial_derivative_over_radius
      use int_c,            only: tf ,tg ,xu, xv
      use orb_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: j
      integer , intent(in) :: jp
      real(real64) , intent(inout) :: p(max_radial_points)
      real(real64) , intent(inout) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: dkhb2, dkhb2f, cpi, cmi, pi, qi, tfi, tgi, cpim1, cmim1, &
         pim1, qim1, tfim1, tgim1, ucim1, udim1, uei
!-----------------------------------------------
!
!   One global initialization
!
      dkhb2 = 0.5_real64*transformed_grid_step*dble(orbital_kappa(j))
!
!   Tabulate P(r) and Q(r) by step-by-step integration
!
!   Initializations: set quantities for I = 3
!
      i = 3
      dkhb2f = dkhb2*radial_derivative_over_radius(i)
      cpi = 1.0_real64 + dkhb2f
      cmi = 1.0_real64 - dkhb2f
      pi = p(i)
      qi = q(i)
      tfi = tf(i)
      tgi = tg(i)
!
!   March out to from I = 4 to I = JP
!
!XHH Use doo-loop
      do i = 4, jp
         cpim1 = cpi
         cmim1 = cmi
         dkhb2f = dkhb2*radial_derivative_over_radius(i)
         cpi = 1.0_real64 + dkhb2f
         cmi = 1.0_real64 - dkhb2f
         pim1 = pi
         qim1 = qi
         tfim1 = tfi
         tgim1 = tgi
         tfi = tf(i)
         tgi = tg(i)
         ucim1 = cmim1*pim1 - tfim1*qim1 + xu(i-1)
         udim1 = cpim1*qim1 - tgim1*pim1 + xv(i-1)
         uei = cpi*cmi - tfi*tgi
         pi = (cmi*ucim1 - tfi*udim1)/uei
         qi = (cpi*udim1 - tgi*ucim1)/uei
         p(i) = pi
         q(i) = qi
      end do
!
      return
      end subroutine out

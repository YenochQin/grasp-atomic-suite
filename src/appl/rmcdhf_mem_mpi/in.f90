!***********************************************************************
!                                                                      *
      subroutine in(iorb, jp, p, q, mtp)
!                                                                      *
!   This program computes the solution of an  inhomogeneous  pair of   *
!   radial Dirac equations in the tail region. A simple extension of   *
!   the method of  C Froese Fischer [C Froese, Can J Phys, 41 (1963)   *
!   1895]  is  used. The  equations  are treated as a boundary value   *
!   problem, with  the value of P(r) given at the inner boundary and   *
!   required to  be  sufficiently small  at the outer boundary.  The   *
!   location  of this second boundary is determined in the course of   *
!   the calculation. The same  finite  difference equations are used   *
!   as in the outward integration.                                     *
!                                                                      *
!   Arguments:                                                         *
!                                                                      *
!      IORB: (Input) Index of orbital                                  *
!      JP:   (Input) Tabulation point at which the outward             *
!            integration terminated                                    *
!      P,Q:  (Input and Output) arrays containing, respectively, the   *
!            large and small components of the solution                *
!      MTP:  (Output) Maximum tabulation point of functions P, Q       *
!                                                                      *
!   When  written  in  matrix form the system of linear equations iS   *
!   M*W = V, where  W  is a vector consisting of the elements Q(JP),   *
!   P(JP+1), Q(JP+1), ... . The matrix  M  is  of band type, with  5   *
!   elements  in each row.  M is  expressed as a product of two tri-   *
!   angular matrices L and  U (the L/U decomposition, carried out by   *
!   the Crout procedure), variable elements of L being stored in the   *
!   arrays TC, TD  and TE, and  variable elements of U in TH, TI and   *
!   TJ. The  solution vector W  is then obtained by solving two tri-   *
!   angular systems, L*Z = V AND U*W = Z.                              *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT.                                       *
!                                           Last update: 10 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: max_radial_points
      use def_c,            only: accy
      use grid_c,           only: transformed_grid_step, radial_point_count, radial_grid, radial_derivative_over_radius
      use int_c,            only: tf, tg, xu, xv
      use orb_c,            only: orbital_kappa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_i, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iorb
      integer , intent(in) :: jp
      integer , intent(out) :: mtp
      real(real64) , intent(inout) :: p(max_radial_points)
      real(real64) , intent(inout) :: q(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, lcnum, mtpp1
      real(real64) :: th(max_radial_points), ti(max_radial_points), tj(max_radial_points), xr(max_radial_points), &
         xs(max_radial_points)
      real(real64) :: eps, hhk, hhkf, cpj, cmj, cpjp1, cmjp1, tei, tthis, tci, &
         tdi, tlast
      character :: cnum*5
!-----------------------------------------------
!
!   Global initializations
!
!ww      EPS = 0.1D 00*ACCY
      eps = 0.01_real64*accy
      hhk = 0.5_real64*transformed_grid_step*dble(orbital_kappa(iorb))
!
!   Initialize counters
!
      i = 1
      j = jp
!
!   Other initializations
!
      hhkf = hhk*radial_derivative_over_radius(j)
      cpj = 1.0_real64 + hhkf
      cmj = 1.0_real64 - hhkf
      hhkf = hhk*radial_derivative_over_radius(j+1)
      cpjp1 = 1.0_real64 + hhkf
      cmjp1 = 1.0_real64 - hhkf
!
!   Compute required elements of first two rows of  L  and  U
!
      th(i) = cpj
      tei = -tf(j)/th(i)
      ti(i) = (-cpjp1) + tei*tg(j+1)
      tj(i) = cmjp1*tei - tf(j+1)
!
!   First elements of solution vector Z
!
      xr(i) = (-xv(j)) + tg(j)*p(j)
      xs(i) = (-xu(j)) - cmj*p(j) - tei*xr(i)
!
      tthis = abs(xs(i)/ti(i))
!
    1 continue
      i = i + 1
      j = j + 1
!
!   Failure if tables not long enough
!
      if (j >= radial_point_count) then
         call convrt (radial_point_count, cnum, lcnum)
         write (6, *) "IN: maximum tabulation point exceeds"
         write (6, *) " dimensional limit (currently "//cnum(1:lcnum)//");"
         write (6, *) " radial wavefunction may indicate a"
         write (6, *) " continuum state."
         stop
      end if
!
!   Compute required elements of remaining rows of  L  and  U
!
      cpj = cpjp1
      cmj = cmjp1
      hhkf = hhk*radial_derivative_over_radius(j+1)
      cpjp1 = 1.0_real64 + hhkf
      cmjp1 = 1.0_real64 - hhkf
!
      tci = -tg(j)/ti(i-1)
      th(i) = cpj - tci*tj(i-1)
      tdi = cmj/ti(i-1)
      tei = ((-tf(j))-tdi*tj(i-1))/th(i)
      ti(i) = (-cpjp1) + tei*tg(j+1)
      tj(i) = (-tf(j+1)) + cmjp1*tei
!
!   Solution of  L*Z = V
!
      xr(i) = (-xv(j)) - tci*xs(i-1)
      xs(i) = (-xu(j)) - tdi*xs(i-1) - tei*xr(i)
!
!   Test for outer boundary
!
      tlast = tthis
      tthis = abs(xs(i)/ti(i))
      if (tthis + tlast <= eps) then
         mtp = j
      else
         goto 1
      end if
!
!   Reset counter
!
      i = i - 1
!
!   Last two rows of solution of  U*W = Z ; evaluation of Q(J)
!
      q(j) = 0.0_real64
      p(j) = xs(i)/ti(i)
      q(j-1) = (xr(i)+tg(j)*p(j))/th(i)
!
!   Solution of  U*W = Z
!
    2 continue
      j = j - 1
      i = i - 1
!
      if (i > 0) then
         p(j) = (xs(i)-tj(i)*q(j))/ti(i)
         q(j-1) = (xr(i)+(1.0_real64-hhk*radial_derivative_over_radius(j))*q(j)+tg(j)*p(j))/th(i)
         goto 2
      end if
!
!   Complete tables with zeroes
!
      mtpp1 = mtp + 1
      p(mtpp1:radial_point_count) = 0.0_real64
      q(mtpp1:radial_point_count) = 0.0_real64
!
      return
      end subroutine in

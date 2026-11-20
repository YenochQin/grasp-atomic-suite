!***********************************************************************
!                                                                      *
      subroutine in(orbital_index, matching_grid_index, large_component_values, small_component_values, mtp)
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
      use grasp_parameters,    only: max_radial_points
      use atomic_calculation_state,            only: accy
      use radial_grid_state,           only: transformed_grid_step, radial_point_count, radial_grid, radial_derivative_over_radius
      use radial_dirac_workspace,            only: tf, tg, xu, xv
      use orbital_csf_basis,            only: orbital_kappa
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: orbital_index
      integer , intent(in) :: matching_grid_index
      integer , intent(out) :: mtp
      real(real64) , intent(inout) :: large_component_values(max_radial_points)
      real(real64) , intent(inout) :: small_component_values(max_radial_points)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, text_length, mtpp1
      real(real64) :: th(max_radial_points), ti(max_radial_points), tj(max_radial_points), xr(max_radial_points), &
         xs(max_radial_points)
      real(real64) :: eps, hhk, hhkf, cpj, cmj, cpjp1, cmjp1, tei, tthis, tci, &
         tdi, tlast
      character :: number_text*5
!-----------------------------------------------
!
!   Global initializations
!
!ww      EPS = 0.1D 00*ACCY
      eps = 0.01_real64*accy
      hhk = 0.5_real64*transformed_grid_step*dble(orbital_kappa(orbital_index))
!
!   Initialize counters
!
      i = 1
      j = matching_grid_index
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
      xr(i) = (-xv(j)) + tg(j)*large_component_values(j)
      xs(i) = (-xu(j)) - cmj*large_component_values(j) - tei*xr(i)
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
         call convrt (radial_point_count, number_text, text_length)
         write (6, *) "IN: maximum tabulation point exceeds"
         write (6, *) " dimensional limit (currently "//number_text(1:text_length)//");"
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
      small_component_values(j) = 0.0_real64
      large_component_values(j) = xs(i)/ti(i)
      small_component_values(j-1) = (xr(i)+tg(j)*large_component_values(j))/th(i)
!
!   Solution of  U*W = Z
!
    2 continue
      j = j - 1
      i = i - 1
!
      if (i > 0) then
         large_component_values(j) = (xs(i)-tj(i)*small_component_values(j))/ti(i)
         small_component_values(j-1) = (xr(i)+(1.0_real64-hhk*radial_derivative_over_radius(j))*small_component_values(j)+tg(j) &
            & *large_component_values(j))/th(i)
         goto 2
      end if
!
!   Complete tables with zeroes
!
      mtpp1 = mtp + 1
      large_component_values(mtpp1:radial_point_count) = 0.0_real64
      small_component_values(mtpp1:radial_point_count) = 0.0_real64
!
      return
      end subroutine in

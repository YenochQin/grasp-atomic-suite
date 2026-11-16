!***********************************************************************
!                                                                      *
      subroutine yzk(k, i, j)
!                                                                      *
!   This subroutine evaluates Hartree Y- and Z-functions:              *
!                                                                      *
!               (K)            (K)           (K)                       *
!              Y   (I,J;r) =  Z   (I,J;r) + W   (I,J;r)                *
!                                                                      *
!   where                                                              *
!                                                                      *
!    (K)                                                               *
!   Z   (I,J;r) =  I ( (s/r)   (P (s)*P (s) + Q (s)*Q (s)) ; 0 - r )   *
!                                I     J       I     J                 *
!                                                                      *
!   and                                                                *
!                                                                      *
!    (K)                    K+1                                        *
!   W   (I,J;r) =  I ( (r/s)   (P (s)*P (s) + Q (s)*Q (s)) ; r -       *
!                                I     J       I     J    INFINITY )   *
!                                                                      *
!   where  I ( G(r,s) ; range )  denotes the integral of G(r,s) over   *
!   range  in  s .  The Y-function is tabulated in  COMMON/TATB/  in   *
!   array  TB , the Z-function in array tA .                           *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!   Modified by Anders Ynnerman, at Vanderbilt         : 03 Feb 1994   *
!   Modified by Jacek Bieron     at Vanderbilt         : 08 Feb 1994   *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:01   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: radial_workspace_capacity
      use cnc_c,           only: cnc5c
      use def_c,           only: accy
      use grid_c
      use ncc_c
      use orb_c
      use tatb_c,          only: zk=>ta, yk=>tb, mtp
      use wave_c,          only: orbital_last_grid_index , radial_large_components, radial_small_components
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: k
      integer  :: i
      integer  :: j
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: klimit = 20
      integer, parameter :: klimit1 = klimit + 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, mtpp1, mtpp3, mtpp4, kk, np4
      real(double), dimension(radial_workspace_capacity) :: rhop
      real(double), dimension(radial_workspace_capacity,klimit) :: rttk, rttkm, rttk1, rttkm1
      real(double), dimension(radial_workspace_capacity) :: rm, wk, temp
      real(double) :: sum, rtmp, zklim, dif
      logical, dimension(0:klimit) :: kcalc

      save kcalc, rttk, rttkm, rttk1, rttkm1, rm
!-----------------------------------------------
!
!      P_OINTER (PNTRPF,PF(NNNP,1)),(PNTRQF,QF(NNNP,1))
!
!
      data kcalc/ klimit1*.false./
      !EQUIVALENCE (TA(1), ZK(1)), (TB(1), YK(1))
      if (.not.kcalc(k)) then
         kcalc(k) = .true.
         if (k > klimit) then
            write (6, *) ' increase klimit in yzk.f '
            stop
         endif
         if (k > 0) then
            rttk(2:radial_point_count,k) = radial_grid(2:radial_point_count)**k
            rttkm(2:radial_point_count,k) = 1.d0/rttk(2:radial_point_count,k)
            rttk1(2:radial_point_count,k) = rttk(2:radial_point_count,k)*radial_grid(2:radial_point_count)
            rttkm1(2:radial_point_count,k) = 1.d0/rttk1(2:radial_point_count,k)
            rm(2:radial_point_count) = 1.d0/radial_grid(2:radial_point_count)
         else
            rm(2:radial_point_count) = 1.d0/radial_grid(2:radial_point_count)
         endif
      endif
!
!   Determine maximum tabulation point as location beyond which
!   RHOP  (see comment statements below) would be zero; determine
!   other important locations
!
      mtp = min(orbital_last_grid_index(i),orbital_last_grid_index(j))
      mtpp1 = mtp + 1
      mtpp3 = mtp + 3
      mtpp4 = mtp + 4
!
!   Compute RP(s)*(P (s)*P (s)+Q (s)*Q (s)) and store in RHOP
!                   I     J     I     J
!
      do ii = 2, mtp
         rhop(ii) = radial_grid_derivative(ii)*(radial_large_components(ii,i)*radial_large_components(ii,j) + &
            & radial_small_components(ii,i)*radial_small_components(ii,j))
      end do
!
!   Fill array TEMP with R**K * RHOP
!
      temp(1) = 0.0d00
      if (k == 0) then
         temp(2:mtp) = rhop(2:mtp)
      else
         temp(2:mtp) = rttk(2:mtp,k)*rhop(2:mtp)
      endif
!
!   Set an additional four points to zero
!
      temp(mtpp1:mtpp4) = 0.0d00
!
!                                     K
!   Compute the first few values of  R  * ZK  using semi-open
!   Newton-Cotes formulae
!
      zk(1) = 0.0d00
      do ii = 2, 4
         sum = 0.0d00
         do kk = 2, 5
            sum = sum + cnc5c(kk,ii)*temp(kk)
         end do
         zk(ii) = sum
      end do
!                         K
!   Compute remainder of R  * ZK: march out to MTP+3; use closed
!   Newton-Cotes formula
!
      do ii = 5, mtpp3
         rtmp = c1*(temp(ii-4)+temp(ii)) + c2*(temp(ii-3)+temp(ii-1)) + c3*temp&
            (ii-2)
         zk(ii) = zk(ii-4) + rtmp
      end do
!
!                                       K   (K)
!   Determine the asymptotic value of  R * Z
!
!                   (0)
!   Correction to  z   : in the manner of  C Froese Fischer,
!   The Hartree-Fock Method for Atoms, John Wiley & Sons,
!   New York, 1977, p 235.
!
      if (k == 0) then
!
         if (i == j) then
            zklim = 1.0d00
         else
            zklim = 0.0d00
         endif
!
         do kk = mtpp3, mtp, -1
            dif = zk(kk) - zklim
            if (abs(dif) <= accy) cycle
            zk(kk:kk-((-2-kk)/(-4)-1)*4:(-4)) = zk(kk:kk-((-2-kk)/(-4)-1)*4:&
               (-4)) - dif
         end do
!
      else
!
         zklim = zk(mtpp3)
!
      endif
!
!   Tabulate  ZK  for entire internal grid
!
      if (k == 0) then
!
         zk(mtpp4:radial_point_count) = zklim
!
      else
!
         zk(2:mtpp3) = zk(2:mtpp3)*rttkm(2:mtpp3,k)
!
         zk(mtpp4:radial_point_count) = zklim*rttkm(mtpp4:radial_point_count,k)
!
      endif
!
!   Start array WK / R**(K+1)
!
      np4 = radial_point_count + 4
      wk(np4:mtpp1:(-1)) = 0.0d00
!
!   Fill array TEMP with RHOP / R**(K+1) ; set TEMP(1) = 0
!   to avoid 0/0 case
!
      temp(1) = 0.0d00
      if (k == 0) then
         temp(2:mtp) = rhop(2:mtp)*rm(2:mtp)
      else
         temp(2:mtp) = rhop(2:mtp)*rttkm1(2:mtp,k)
      endif
!
!   Compute remainder of WK / R**(K+1): march in to the origin
!
      do ii = mtp, 2, -1
         wk(ii) = wk(ii+4) + c1*(temp(ii)+temp(ii+4)) + c2*(temp(ii+1)+temp(ii+&
            3)) + c3*temp(ii+2)
      end do
      wk(1) = 0.0d00
!
!   Compute WK
!
      if (k == 0) then
         wk(2:mtp) = wk(2:mtp)*radial_grid(2:mtp)
      else
         wk(2:mtp) = wk(2:mtp)*rttk1(2:mtp,k)
      endif
!
!   Assemble solution
!
      yk(1) = 0.0d00
      yk(2:radial_point_count) = zk(2:radial_point_count) + wk(2:radial_point_count)
!
      return
      end subroutine yzk

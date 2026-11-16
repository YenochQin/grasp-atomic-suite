!***********************************************************************
!                                                                      *
      subroutine zkf(k, i, j)
!                                                                      *
!   This subroutine evaluates Hartree Z-functionals:                   *
!                                                                      *
!              (k)                     k                               *
!             Z   [f(r);r] =  I ( (s/r)   f(s) ; 0 - r )               *
!                                                                      *
!   where  I ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
!   range  in  s .    The Z-functional is tabulated in  COMMON/TATB/   *
!   in array  TB . The f-function is assumed tabulted in array  TA .   *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 14 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: radial_workspace_capacity
      use cnc_c, only: cnc5c
      use grid_c, only: radial_point_count, radial_grid, radial_grid_derivative
      use ncc_c
      use tatb_c, zk=>tb
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: k
      integer  :: i
      integer  :: j    !!!! Arument not referenced !!!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, mtpp1, mtpp3, mtpp4, kk
      real(double), dimension(radial_workspace_capacity) :: rhop, rttk, temp
      real(double) :: sum, zklim
!-----------------------------------------------
!
      if (k > 0) then
         rttk(2:radial_point_count) = radial_grid(2:radial_point_count)**k
      endif
!
!   MTP is fed in through  COMMON/TATB/
!
      mtpp1 = mtp + 1
      mtpp3 = mtp + 3
      mtpp4 = mtp + 4
!
!   Compute  RP(S)*F(S)  and store it in  RHOP
!
      rhop(2:mtp) = radial_grid_derivative(2:mtp)*ta(2:mtp)
!
!   Fill array TEMP with r**k * RHOP
!
      temp(1) = 0.0d00
      if (k == 0) then
         temp(2:mtp) = rhop(2:mtp)
      else
         temp(2:mtp) = rttk(2:mtp)*rhop(2:mtp)
      endif
!
!   Set an additional four points to zero
!
      temp(mtpp1:mtpp4) = 0.0d00
!
!                                     k
!   Compute the first few values of  r  * ZK  using semi-open
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
!                         k
!   Compute remainder of r  * ZK: march out to MTP+3
!
      do ii = 5, mtpp3
         zk(ii) = zk(ii-4) + c1*(temp(ii-4)+temp(ii)) + c2*(temp(ii-3)+temp(ii-&
            1)) + c3*temp(ii-2)
      end do
!                                       k   (k)
!   Determine the asymptotic value of  r * Z
!
!   Compute ZK
!
      zklim = zk(mtpp3)
!
      if (k == 0) then
!
         zk(mtpp4:radial_point_count) = zklim
!
      else
!
         zk(2:mtpp3) = zk(2:mtpp3)/rttk(2:mtpp3)
!
         zk(mtpp4:radial_point_count) = zklim/rttk(mtpp4:radial_point_count)
!
      endif
!
      return
      end subroutine zkf

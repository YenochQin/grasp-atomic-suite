!***********************************************************************
!                                                                      *
      subroutine quad(result)
!                                                                      *
!   The argument result is an approximation  to the integral of F(R)   *
!   from  zero  to  infinity,  where the values of RP(I)*F(R(I)) are   *
!   tabulated in the array  TA(I). The integral in the interval zero   *
!   to R(J) is computed by use of an analytical fit                    *
!                                                                      *
!                                SIGMA                                 *
!                      F(R) = A R                                      *
!                                                                      *
!   A five-point  Closed  Newton-Cotes  formula (cf. F B Hildebrand,   *
!   Introduction to Numerical Analysis, second edition, McGraw-Hill,   *
!   New York, 1974, p 93)  is  used  to  compute the integral in the   *
!   interval  R(J:MTP).  The  contribution  from  the  tail  of  the   *
!   function beyond the last  tabular  point  (MTP) is assumed to be   *
!   negligible. The method uses  MTP+3  tabulation points. Array  TA   *
!   should therefore be dimensioned to at least  N+4 .                 *
!                                                                      *
!   No subroutines called.                                             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:18   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use def_c
      use grid_c
      use ncc_c
      use tatb_c, only: ta, mtp
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double) , intent(out) :: result
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mtpm1, i, ip1, loc
      real(double) :: tai, taip1, quott, frip1, fri, ratio, rip1, ri, sigma
!-----------------------------------------------
!
!
!   Find first values that will permit computation of exponent
!
      mtpm1 = mtp - 1
      do i = 2, mtpm1
!
         tai = ta(i)
         if (abs(tai) <= 0.d0) cycle
!
         ip1 = i + 1
         taip1 = ta(ip1)
         quott = taip1/tai
!
         if (quott <= 0.d0) cycle
!
!   Exponent from fit
!
         frip1 = taip1/rp(ip1)
         fri = tai/rp(i)
         ratio = frip1/fri
         rip1 = r(ip1)
         ri = r(i)
         sigma = log(ratio)/log(rip1/ri)
!
!   Analytical integration and error estimate for interval r(1:i)
!
         fri = ri*fri
         result = fri/(sigma + 1.d0)
!
!   Set the tail to zero
!
         ta(mtp+1:3+mtp) = 0.d0
!
!   Newton-Cotes quadature for the remainder
!
         result = result + c1*tai
         result = result + sum(c2*(ta(ip1:mtp:4)+ta(ip1+2:mtp+2:4))+c3*ta(ip1+1&
            :mtp+1:4)+c4*ta(ip1+3:mtp+3:4))
         if (mod(mtp - i,4) == 0) result = result - c1*ta(mtp)
!
!   Test of result's accuracy; `decomment' to activate
!
!              ESTDER = 10.0D 00*RI*FRI
!              RATIO = ABS (ESTDER/RESULT)
!              IF (RATIO .GT. ACCY) PRINT (*,300) RATIO
!
         go to 4
!
      end do
!
!   No value which will permit computation of exponent
!
      result = 0.d0
!
    4 continue
      return
!
  300 format(/,'QUAD: Estimated accuracy is ',1p,d10.3,/,&
         ' Decrease RNT or improve input data conditioning to',' ameliorate.'/)
      return
!
      end subroutine quad

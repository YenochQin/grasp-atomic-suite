!***********************************************************************
!                                                                      *
      real(real64) function funl (x, k)
!                                                                      *
!   This  function  evaluates the LK(X) functions using the analytic   *
!   functions defined  in table 5  and equations  (20) and  (21)  of   *
!   Fullerton and Rinker.                                              *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: k
      real(real64), intent(in) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k1, i
      real(real64) :: f(6, 2)
      real(real64) :: g(3, 2)
      real(real64) :: h(2, 2)
      real(real64) :: a, b, sum, xn, x2, sumg, sumh, xm
!-----------------------------------------------
!
!
      data f/ 2.008188_real64, -2.397605_real64, 1.046471_real64, -3.670660e-1_real64, &
         6.374000e-2_real64, -3.705800e-2_real64, 1.646407_real64, -2.092942_real64, 9.623100e-1_real64, &
         -2.549600e-1_real64, 1.644040e-1_real64, 0.0_real64/
!
      data g/ 7.51198e-1_real64, 1.38889e-1_real64, 2.0886e-2_real64, 1.37691e-1_real64, -4.16667e-1_real64, &
          - 9.7486e-2_real64/
!
      data h/ -4.44444e-1_real64, -3.472e-3_real64, 4.44444e-1_real64, 1.7361e-2_real64/
!
      data a, b/ 2.2_real64,  - 1.72_real64/
!
      if (k>=0 .and. k<=1) then
         if (x <= 2.0_real64) then
            if (x == 0.0_real64) goto 6
!
!   Use rational approximation for X < 2
!
            k1 = k + 1
            sum = 0.0_real64
            xn = 1.0_real64
            do i = 1, 6
               sum = sum + xn*f(i,k1)
               xn = xn*x
            end do
            x2 = x*x
            sumg = g(1,k1) + x2*(g(2,k1)+x2*g(3,k1))
            sumh = h(1,k1) + x2*x2*h(2,k1)
            xn = log(x)
            sumg = xn*(sumg + xn*sumh)
            if (k /= 0) then
               sum = sum + sumg
               goto 7
            end if
            sum = sum + x*sumg
            goto 7
         end if
!
         sum = a + b/x
         if (k /= 0) sum = sum + (sum + b/x)/x
         sum = sum/x
         xm = -x
         sum = sum*exp(xm)
         goto 7
    6    continue
         if (k == 1) goto 98
         sum = f(1,1)
    7    continue
         funl = sum
         return
!
!   Error section
!
   98    continue
         write (*, 302)
         stop
      end if
      write (*, 301)
      stop
!
  301 format(/,"FUNL: K must be either 0 or 1")
  302 format(/,"FUNL: Attempt to calculate function for"/,&
         " zero argument and K value of 1")
      return
!
      end function funl

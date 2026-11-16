!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function funl (x, k)
!                                                                      *
!   This  function  evaluates the LK(X) functions using the analytic   *
!   functions defined  in table 5  and equations  (20) and  (21)  of   *
!   Fullerton and Rinker.                                              *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: k
      real(double), intent(in) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k1, i
      real(double), dimension(6,2) :: f
      real(double), dimension(3,2) :: g
      real(double), dimension(2,2) :: h
      real(double) :: a, b, sum, xn, x2, sumg, sumh, xm
!-----------------------------------------------
!
!
      data f/ 2.008188d00, -2.397605d00, 1.046471d00, -3.670660d-01, &
         6.374000d-02, -3.705800d-02, 1.646407d00, -2.092942d00, 9.623100d-01, &
         -2.549600d-01, 1.644040d-01, 0.0d00/
!
      data g/ 7.51198d-01, 1.38889d-01, 2.0886d-02, 1.37691d-01, -4.16667d-01, &
          - 9.7486d-02/
!
      data h/ -4.44444d-01, -3.472d-03, 4.44444d-01, 1.7361d-02/
!
      data a, b/ 2.2d00,  - 1.72d00/
!
      if (k>=0 .and. k<=1) then
         if (x <= 2.0d00) then
            if (x == 0.0d00) go to 6
!
!   Use rational approximation for X < 2
!
            k1 = k + 1
            sum = 0.0d00
            xn = 1.0d00
            do i = 1, 6
               sum = sum + xn*f(i,k1)
               xn = xn*x
            end do
            x2 = x*x
            sumg = g(1,k1) + x2*(g(2,k1)+x2*g(3,k1))
            sumh = h(1,k1) + x2*x2*h(2,k1)
            xn = dlog(x)
            sumg = xn*(sumg + xn*sumh)
            if (k /= 0) then
               sum = sum + sumg
               go to 7
            endif
            sum = sum + x*sumg
            go to 7
         endif
!
         sum = a + b/x
         if (k /= 0) sum = sum + (sum + b/x)/x
         sum = sum/x
         xm = -x
         sum = sum*dexp(xm)
         go to 7
    6    continue
         if (k == 1) go to 98
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
      endif
      write (*, 301)
      stop
!
  301 format(/,'FUNL: K must be either 0 or 1')
  302 format(/,'FUNL: Attempt to calculate function for'/,&
         ' zero argument and K value of 1')
      return
!
      end function funl

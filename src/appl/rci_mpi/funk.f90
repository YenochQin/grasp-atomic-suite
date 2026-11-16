!***********************************************************************
!                                                                      *
      real(kind(0.0d0)) function funk (x, n)
!                                                                      *
!   This  function  evaluates the KN(X) functions using the analytic   *
!   functions defined in tables 1 and 3 of Fullerton and Rinker.       *
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
      integer, intent(in) :: n
      real(double), intent(in) :: x
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer, dimension(4) :: np
      integer :: k, nn, i
      real(double), dimension(10,4) :: p
      real(double), dimension(2,4) :: b
      real(double), dimension(3,4) :: c
      real(double), dimension(5,4) :: d, e
      real(double) :: xn, sum, x2, bsum, csum, dsum, esum, xm
!-----------------------------------------------
!
!
      data p/ 8.8357293375d-01, -2.8259817381d-01, -5.8904879578d-01, &
         1.2500133434d-01, -3.2729913852d-02, 8.2888574511d-03, &
         -1.0327765800d-05, 6.3643668900d-05, 0.0d00, 0.0d00, -7.1740181754d-01&
         , 1.1780972274d00, -3.7499963087d-01, 1.3089675530d-01, &
         -3.8258286439d-02, -2.4297287300d-05, -3.5920148670d-04, &
         -1.7170090700d-05, 0.0d00, 0.0d00, 9.9999999987d-01, 1.9770200000d-08&
         , -7.5000050190d-01, 7.8540306316d-01, -3.4988601655d-01, &
         6.4596333000d-05, -9.8189080747d-03, 8.6513145800d-05, &
         -2.3969236620d-04, 0.0d00, 6.0000000002d00, -6.4305200000d-08, &
         2.1049413000d-06, -2.6711271500d-05, -1.3705236152d-01, &
         -6.3476104090d-04, -7.8739801501d-02, -1.9641740173d-03, &
         -3.4752369349d-03,  - 7.3145316220d-04/
!
      data b/ -3.19999594323d+02, 2.53900995981d00, -6.40514843293d+01, &
         7.11722714285d-01, 5.19010136460d+03, 8.28495496200d+01, &
         3.18150793824d+02, 4.33898867347d+01/
!
      data c/ -3.19999594333d+02, 2.53901020662d+00, 0.0d00, 6.40514843287d+01&
         , -7.11722686403d-01, 8.042207748d-04, 2.76805406060d+04, &
         -3.27039477790d+02, 0.0d00, 8.48402116837d+02, -2.56939867765d+01, &
         3.20844906346d-01/
!
      data d/ 5.018065179d+00, 7.151891262d+01, 2.116209929d+02, &
         3.140327478d+01, -1.0d00, 2.172386409d+02, 1.643364528d+03, &
         2.122244512d+03, -4.512004044d+01, 1.0d00, 8.540770444d00, &
         6.076242766d+01, 9.714630584d+01, 3.154973593d+01, 1.0d00, &
         5.9243015865d-01, 2.0596312871d00, 3.7785190424d00, 3.5614853214d00, &
         1.0d00/
!
      data e/ 2.669207401d+00, 5.172549669d+01, 2.969809720d+02, &
         5.364324164d+02, 1.535335924d+02, 1.155589983d+02, 1.292191441d+03, &
         3.831198012d+03, 2.904410075d+03, 0.0d00, 4.543392478d00, &
         3.514920169d+01, 6.019668656d+01, 8.468839579d00, 0.0d00, &
         3.1511867816d-01, 3.473245222d-01, 3.8791936870d-02, -1.3059741497d-03&
         , 0.0d00/
!
      data np/ 8, 8, 9, 10/
!
      if (x /= 0.0d00) then
         if (n<0 .or. n==2 .or. n==4 .or. n>5) go to 98
         if (n - 3 <= 0) then
            if (n - 3 /= 0) then
               k = n + 1
               xn = 1.0d00
               go to 4
            endif
            k = n
            xn = 1.0d00/x**2
            go to 4
         endif
         k = n - 1
         xn = 1.0d00/x**4
    4    continue
         if (x <= 1.0d00) then
!
!   Calculate function for X < = 1
!
            nn = np(k)
            sum = 0.0d00
            do i = 1, nn
               sum = sum + p(i,k)*xn
               xn = xn*x
            end do
            x2 = x*x
            bsum = b(1,k) + x2*(b(2,k)+x2)
            csum = c(1,k) + x2*(c(2,k)+x2*c(3,k))
            go to (6,8,7,8) k
    6       continue
            bsum = bsum*x
            go to 8
    7       continue
            bsum = bsum*x2
    8       continue
            sum = sum + bsum*dlog(x)/csum
            funk = sum
            return
         endif
!
!   Calculate function for X > 1
!
         xn = 1.0d00
         dsum = 0.0d00
         esum = 0.0d00
         do i = 1, 5
            dsum = dsum + d(i,k)*xn
            esum = esum + e(i,k)*xn
            xn = xn/x
         end do
         xm = -x
         sum = dsum*exp(xm)/(esum*dsqrt(x**3))
         funk = sum
         return
      endif
      if (n /= 0) go to 99
      funk = p(1,1)
      return
!
!   Error section
!
   98 continue
      write (*, 302)
      stop
   99 continue
      write (*, 301)
      stop
!
  301 format(/,' Attempt to calculate FUNK (0,N) for N > 0')
  302 format(/,' Attempt to calculate FUNK (X,N) for N other than',&
         ' 0, 1, 3 and 5.')
      return
!
      end function funk

!***********************************************************************
!                                                                      *
      subroutine cxk(s, is, kaps, nu, k, ibr, iex)
!                                                                      *
!   Computes  the  coefficients of radial integrals in the expansion   *
!   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS  : Orbital labels                                            *
!      KAPS: Values of 2*kappa                                         *
!      NU  : Order of radial integral                                  *
!      K   : Index of tensor operator                                  *
!      IEX : 1 for direct, 2 for exchange terms                        *
!      IBR : Classifies type of radial integral.There are 4 distinct   *
!            cases:                                                    *
!            IBR = 1 A. All states distinct                            *
!                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
!                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
!                    These give 12 distinct  radial  integrals, with   *
!                    values of K and NU limited only by angular mom-   *
!                    entum and parity                                  *
!            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
!                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
!                    This case gives one non-zero integral when K =    *
!                    NU is ODD                                         *
!            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
!                    (IA .NE. IB)                                      *
!                    Integrals of magnetic F-type when K = NU is odd   *
!            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
!                    netic G-type integrals and  four  H-TYPE  inte-   *
!                    grals                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!      S   : Coefficients S(MU) MU = 1,12                              *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92] CRE.                                           *
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
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cre_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nu
      integer, intent(in) :: k
      integer, intent(in) :: ibr
      integer, intent(in) :: iex
      integer, intent(in) :: is(4)
      integer, intent(in) :: kaps(4)
      real(real64), intent(inout) :: s(12)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mu, ia, ib, ic, id, ka, kb, kc, kd, kk, ik, ip
      real(real64) :: d, h, dk1, dk2, fk, gk, g1, g2, g3, g4, a, f1, f2, f3, f4&
         , b, dk
!-----------------------------------------------
!
!
!   1.0  Initialization
!
      s = 0.0_real64
!
      ia = is(1)
      ib = is(2)
      ic = is(3)
      id = is(4)
      ka = kaps(1)/2
      kb = kaps(2)/2
      kc = kaps(3)/2
      kd = kaps(4)/2
      if (iex == 2) then
         kk = kd
         ik = id
         kd = kc
         id = ic
         kc = kk
         ic = ik
      endif
      select case (ibr)
      case default
         go to 17
!
!   2.0  IBR = 1 --- The general case
!
      case (1)
         if (nu - k >= 0) then
            if (nu - k <= 0) then
               s(1) = -(ka + kc)*(kd + kb)
               if (k == 0) go to 16
               d = k*(k + 1)
               h = cre(ka,k,kc)*cre(kb,k,kd)
               if (mod(k,2) /= 0) h = -h
               s(1) = s(1)*h/d
               s(2:4) = s(1)
               return
            endif
!
!   2.2  NU = K+1
!
            dk1 = kc - ka
            dk2 = kd - kb
            fk = k
            gk = k + 1
            g1 = dk1 - gk
            g2 = dk1 + gk
            g3 = dk2 - gk
            g4 = dk2 + gk
            kk = k + k + 1
            h = cre(ka,k,kc)*cre(kb,k,kd)
            if (mod(k,2) /= 0) h = -h
            a = h*fk/gk/dble(kk*(kk + 2))
            s(1) = a*g1*g3
            s(2) = a*g2*g4
            s(3) = a*g1*g4
            s(4) = a*g2*g3
            return
         endif
!
!   2.2  NU = K-1
!
         dk1 = kc - ka
         dk2 = kd - kb
         fk = k
         gk = k + 1
         f1 = dk1 - fk
         f2 = dk1 + fk
         f3 = dk2 - fk
         f4 = dk2 + fk
         g1 = dk1 - gk
         g2 = dk1 + gk
         g3 = dk2 - gk
         g4 = dk2 + gk
         kk = k + k + 1
         h = cre(ka,k,kc)*cre(kb,k,kd)
         if (mod(k,2) /= 0) h = -h
         a = h*gk/fk/dble(kk*(kk - 2))
         s(1) = a*f2*f4
         s(2) = a*f1*f3
         s(3) = a*f2*f3
         s(4) = a*f1*f4
         b = h/dble(kk*kk)
         s(5) = b*f2*g3
         s(6) = b*f4*g1
         s(7) = b*f1*g4
         s(8) = b*f3*g2
         s(9) = b*f2*g4
         s(10) = b*f3*g1
         s(11) = b*f1*g3
         s(12) = b*f4*g2
         return
!
!   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
!
      case (2)
         if (ia/=ic .or. ib==id) then
            if (ia==ic .or. ib/=id) go to 17
!
            ik = ib
            ib = ia
            ia = ik
            ik = id
            id = ic
            ic = ik
!
            kk = kb
            kb = ka
            ka = kk
            kk = kd
            kd = kc
            kc = kk
         endif
!
         if (mod(k,2) /= 1) return
         dk = k*(k + 1)
         h = cre(ka,k,kc)*cre(kb,k,kd)/dk
         s(1) = h*dble(4*ka*(kb + kd))
         return
!
!   4.0  IBR = 3. Direct magnetic F-integrals
!
      case (3)
         if (ia/=ic .or. ib/=id) go to 17
         if (mod(k,2) /= 1) return
         dk = k*(k + 1)
         h = cre(ka,k,ka)*cre(kb,k,kb)/dk
         s(1) = h*dble(16*ka*kb)
         return
!
!   5.0   IBR = 4. Exchange magnetic G- and H-integrals
!
      case (4)
         if (ia/=id .or. ib/=ic) go to 17
         if (nu - k >= 0) then
            if (nu - k <= 0) then
               s(1) = dble(ka + kb)*cre(ka,k,kb)
               ip = abs(ka) - abs(kb) + k + 1
               s(1) = s(1)*s(1)/dble(k*(k + 1))
               if (mod(ip,2) /= 0) s(1) = -s(1)
               s(3) = s(1)
               s(2) = s(1) + s(1)
               return
            endif
!
!   5.2  NU = K+1
!
            dk = kb - ka
            gk = k + 1
            fk = k
            g1 = dk + gk
            g2 = dk - gk
            kk = k + k + 1
            h = cre(ka,k,kb)**2
            if (ka*kb < 0) h = -h
            a = h*fk/gk/dble(kk*(kk + 2))
            s(1) = -a*g1*g1
            s(2) = -2.0_real64*a*g1*g2
            s(3) = -a*g2*g2
            return
         endif
!
!   5.3  NU = K-1
!
         dk = kb - ka
         fk = k
         gk = k + 1
         f1 = dk + fk
         f2 = dk - fk
         g1 = dk + gk
         g2 = dk - gk
         kk = k + k + 1
         h = cre(ka,k,kb)**2
         if (ka*kb < 0) h = -h
         a = h*gk/fk/dble(kk*(kk - 2))
         s(1) = -a*f2*f2
         s(2) = -2.0_real64*a*f1*f2
         s(3) = -a*f1*f1
         b = h/dble(kk*kk)
         b = b + b
         s(4) = -b*f1*g2
!     S(5) = S(4)
         s(5) = -b*f2*g1
         s(6) = -b*f1*g1
         s(7) = -b*f2*g2
         return
      end select
!
!   6.0  Special cases and errors
!
!   Illegal zero value of K in Type 1
!
   16 continue
      write (*, 300) is(1), is(2), is(3), is(4), nu, ibr, iex
      stop
!
!   Illegal combination of states in Type 3 or 4
!
   17 continue
      write (*, 301) ibr, is(1), is(2), is(3), is(4), nu, k, iex
      stop
!
  300 format('CXK: Illegal value K = 0 -'/,1x,4i3,2x,i3,2x,2i2)
  301 format('CXK: Type ',i2,'-'/,1x,i2,3x,4i3,2x,2i3,2x,i2)
      return
!
      end subroutine cxk

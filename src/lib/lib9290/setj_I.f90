      module setj_i
      implicit none (type, external)
      interface
      subroutine setj (is, js, ks, ns, kj23)
      implicit none (type, external)
      integer, parameter :: mangm = 60
      integer, parameter :: mtriad = 12
      integer, intent(in) :: is(2, 2)
      integer, intent(in) :: js(2, 2)
      integer, intent(in) :: ks(2, 2)
      integer, intent(in) :: ns
      integer, intent(inout) :: kj23
!VAST.../L1/ JBQ1(IN), JTQ1(IN), JTQ2(IN)
!VAST.../L2/ J2S(INOUT), J3S(INOUT)
!VAST.../M0/ JJC1(IN), JJC2(IN)
!VAST.../M2/ JJQ1(IN), JJQ2(IN)
!VAST.../M3/ JLIST(IN)
!VAST.../COUPLE/ MJA(OUT), NJA(OUT), J1(OUT), J2(INOUT), J3(INOUT)
      end subroutine
      end interface
      end module

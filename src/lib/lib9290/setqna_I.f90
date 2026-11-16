      module setqna_i
      implicit none (type, external)
      interface
      subroutine setqna (ja, jb)
      implicit none (type, external)
      integer, intent(in) :: ja
      integer, intent(in) :: jb
!VAST.../M0/ JJC1(INOUT), JJC2(INOUT)
!VAST.../M1/ NQ1(OUT), NQ2(OUT)
!VAST.../M2/ JJQ1(OUT), JJQ2(OUT)
!VAST.../M3/ JLIST(INOUT), KLIST(OUT), NPEEL(OUT), NCORE(OUT)
!VAST.../ORB2/ NW(IN)
!VAST...Calls: IQ, JQS, ICHOP, JCUP
      end subroutine
      end interface
      end module

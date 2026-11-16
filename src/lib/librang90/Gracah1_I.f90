      module gracah1_i
      interface
!
      subroutine gracah1(i,j,k,l,m,n,rac)
      use vast_kind_param, only:  double
      integer mfact
      parameter (mfact = 500)
      integer , intent(in) :: i
      integer , intent(in) :: j
      integer , intent(in) :: k
      integer , intent(in) :: l
      integer , intent(in) :: m
      integer , intent(in) :: n
      real(double) , intent(out) :: rac
      end subroutine
      end interface
      end module

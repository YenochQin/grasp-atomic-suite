      module ti1tv_i
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine ti1tv(ciin,ncsf,nciv,i,l,t,nshl,ciout,ntestg)
      use vast_kind_param, only:  double
      implicit none
      integer, intent(in) :: ncsf, nciv, nshl, i, l, ntestg
      real(double), dimension(ncsf, nciv), intent(in) :: ciin
      real(double), dimension(nshl), intent(in) :: t
      real(double), dimension(ncsf, nciv), intent(out) :: ciout
      end subroutine
      end interface
      end module

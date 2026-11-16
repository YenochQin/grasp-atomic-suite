      module ti1tv_i
      implicit none (type, external)
      interface
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
      subroutine ti1tv(ciin,ncsf,nciv,i,l,t,nshl,ciout,ntestg)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ncsf, nciv, nshl, i, l, ntestg
      real(real64), intent(in) :: ciin(ncsf, nciv)
      real(real64), intent(in) :: t(nshl)
      real(real64), intent(out) :: ciout(ncsf, nciv)
      end subroutine ti1tv
      end interface
      end module ti1tv_i

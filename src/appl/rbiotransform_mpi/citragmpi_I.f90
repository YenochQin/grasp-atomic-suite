      module citrag_i
      implicit none (type, external)
      interface
      subroutine citrag (ciin,ncsf,nciv,l,nshl,t,nin,ntestg,ciout,scr)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ncsf
      integer, intent(in) :: nciv
      real(real64) :: ciin(ncsf, nciv)
      integer, intent(in) :: l
      integer, intent(in) :: nshl
      real(real64), intent(in) :: t(nshl, nshl)
      integer, intent(in) :: nin
      integer, intent(in) :: ntestg
      real(real64) :: ciout(ncsf, nciv)
      real(real64) :: scr(ncsf, nciv)
      end subroutine citrag
      end interface
      end module citrag_i

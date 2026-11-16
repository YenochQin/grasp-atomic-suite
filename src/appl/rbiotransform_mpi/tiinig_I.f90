      module tiinig_i
      implicit none (type, external)
      interface
      subroutine tiinig (ciin, ncsf, nciv, i, l, const, ciout, ntestg)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: ncsf
      integer, intent(in) :: nciv
      real(real64), intent(in) :: ciin(ncsf, nciv)
      integer, intent(in) :: i
      integer, intent(in) :: l
      real(real64), intent(in) :: const
      real(real64), intent(inout) :: ciout(ncsf, nciv)
      integer, intent(in) :: ntestg
      end subroutine tiinig
      end interface
      end module tiinig_i

      module pamtmt_i
      implicit none (type, external)
      interface
      subroutine pamtmt (x, t, work, norb)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: norb
      real(real64) :: x(norb, norb)
      real(real64), intent(out) :: t(norb, norb)
      real(real64), intent(inout) :: work(*)
      end subroutine pamtmt
      end interface
      end module pamtmt_i

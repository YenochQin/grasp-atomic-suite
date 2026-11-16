      module intrpqi_i
      implicit none (type, external)
      interface
      subroutine intrpqi (pa, qa, ma, ra, j, dnorm)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(in) :: pa(*)
      real(real64), intent(in) :: qa(*)
      real(real64), intent(in) :: ra(*)
      integer, intent(in) :: ma
      integer, intent(in) :: j
      real(real64), intent(out) :: dnorm
      end subroutine
      end interface
      end module

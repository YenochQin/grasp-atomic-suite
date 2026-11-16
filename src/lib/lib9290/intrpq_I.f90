      module intrpq_i
      implicit none (type, external)
      interface
      subroutine intrpq (pa, qa, ma, ra, j, dnorm)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, parameter :: mxord = 13
      real(real64), intent(in) :: pa(*), qa(*), ra(*)
      integer, intent(in) :: ma
      integer, intent(in) :: j
      real(real64), intent(out) :: dnorm
      end subroutine
      end interface
      end module

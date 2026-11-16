      module dampor_i
      implicit none (type, external)
      interface
      subroutine dampor (j, inv, odampj)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(inout) :: inv
      real(real64), intent(in) :: odampj
      end subroutine dampor
      end interface
      end module dampor_i

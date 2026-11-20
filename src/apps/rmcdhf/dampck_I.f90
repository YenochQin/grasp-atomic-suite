      module dampck_i
      implicit none (type, external)
      interface
      subroutine dampck (ipr, orbital_index, ed1, ed2)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(inout) :: ipr
      integer, intent(in) :: orbital_index
      real(real64), intent(inout) :: ed1
      real(real64), intent(inout) :: ed2
      end subroutine dampck
      end interface
      end module dampck_i

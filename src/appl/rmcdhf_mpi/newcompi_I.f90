      module newcompi_i
      implicit none (type, external)
      interface
      subroutine newcompi (sum)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      real(real64), intent(out) :: sum
      end subroutine newcompi
      end interface
      end module newcompi_i

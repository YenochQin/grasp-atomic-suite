      module getaldwt_i
      implicit none (type, external)
      interface
      subroutine getaldwt (csf_count, level_weights)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: csf_count
      real(real64), intent(inout) :: level_weights(csf_count)
      end subroutine getaldwt
      end interface
      end module getaldwt_i

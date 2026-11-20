      module getoldwt_i
      implicit none (type, external)
      interface
      subroutine getoldwt (ndef, ncmin, level_weights)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer :: ndef
!VAST...Dummy argument NDEF is not referenced in this routine.
      integer, intent(in) :: ncmin
      real(real64), intent(inout) :: level_weights(ncmin)
      end subroutine getoldwt
      end interface
      end module getoldwt_i

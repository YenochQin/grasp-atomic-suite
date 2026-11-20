      module engoutgg_i
      implicit none (type, external)
      interface
      subroutine engoutgg(e,level_indices,level_count,mode)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: level_count
      real(real64), intent(in) :: e(level_count)
      integer, intent(in) :: level_indices(level_count)
      integer, intent(in) :: mode
      end subroutine engoutgg
      end interface
      end module engoutgg_i

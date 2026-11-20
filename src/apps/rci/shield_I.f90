      module shield_i
      implicit none (type, external)
      interface
      real(real64) function shield (j)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: j
      end function shield
      end interface
      end module shield_i

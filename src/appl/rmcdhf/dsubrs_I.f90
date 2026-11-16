      module dsubrs_i
      implicit none (type, external)
      interface
      real(real64) function dsubrs (eol, i, j, jblock)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      logical, intent(in) :: eol
      integer, intent(in) :: i
      integer, intent(in) :: j
      integer :: jblock
      end function
      end interface
      end module

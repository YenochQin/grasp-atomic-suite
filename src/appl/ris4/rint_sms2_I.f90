      module rint_sms2_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      real(real64) function rint_sms2 (i, j)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: i
      integer, intent(in) :: j
      end function
      end interface
      end module

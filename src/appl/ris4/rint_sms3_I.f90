      module rint_sms3_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      real(real64) function rint_sms3 (i, j)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: i
      integer, intent(in) :: j
      end function rint_sms3
      end interface
      end module rint_sms3_i

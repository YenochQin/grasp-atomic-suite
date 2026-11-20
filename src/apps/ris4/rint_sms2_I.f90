      module rint_sms2_i
      implicit none (type, external)
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      real(real64) function rint_sms2 (first_orbital, second_orbital)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: first_orbital
      integer, intent(in) :: second_orbital
      end function rint_sms2
      end interface
      end module rint_sms2_i

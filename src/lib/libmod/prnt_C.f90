!
!***********************************************************************
!                                                                      *
      module prnt_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64, byte => int8
      implicit none (type, external)
      private :: real64
      integer :: nvec, nvecmx
      integer :: nvecff, nvecmxff
      integer :: nvecii, nvecmxii
      real(real64) :: pnivecii
      integer, pointer :: ivec(:), ivecff(:), ivecii(:)
      end module prnt_c

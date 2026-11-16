!
!***********************************************************************
!                                                                      *
      module facts_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer, parameter :: mfact = 500
      real(real64) :: gam(mfact)
      end module facts_c

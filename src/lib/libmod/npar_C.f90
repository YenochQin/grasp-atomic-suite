!
!***********************************************************************
!                                                                      *
      module npar_c
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer :: nuclear_model_parameter_count
      real(real64) :: nuclear_model_parameters(2)
      end module npar_c

      module nsmdat_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      real(real64) :: sqn,dmomnm,qmomb
      real(real64) :: hfsi, hfsd, hfsq
      real(real64) :: smsi, smsd, smsq
      end module nsmdat_c

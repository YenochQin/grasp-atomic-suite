!
!***********************************************************************
!                                                                      *
      module numeric_constants
!                                                                      *
!***********************************************************************
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      real(real64) :: zero  = 0.0_real64, &
                      half  = 0.5_real64, &
                      tenth = 0.1_real64, &
                      one   = 1.0_real64, &
                      two   = 2.0_real64, &
                      three = 3.0_real64, &
                      four  = 4.0_real64, &
                      seven = 7.0_real64, &
                      ten   =10.0_real64, &
                      eleven=11.0_real64, &
                      eps   = 1.0e-8_real64
      end module numeric_constants
